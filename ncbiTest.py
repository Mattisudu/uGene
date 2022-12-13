"""
This is a simple script to find ncbi path witch one are effected by the ncbi response architecture, witch one could
be changed or updated with new parameters.
Path means a string with dots, witch one describe a ncbi response object access. Dictionary identifier have to be
separated by dots and list are able to access with square brackets within an index or and index operation.
All in all this path are required for the outToCSV.py -ncbi parameter.
For example here are some id and databases ncbi559295 taxonomy, XP_033769018.1 protein
"""

import pandas as pd
import json
from Bio import Entrez

Entrez.email = "max@muellmail.com"


def getDbInfo(df_name="taxonomy"):
    """
    Gives simply informations about the ncbi database.

    :param str df_name: Name of the ncbi database.
    :return: Returns the info.
    """
    handle = Entrez.einfo(db=df_name)
    record = Entrez.read(handle)
    handle.close()
    return record


def get(id, db="taxonomy"):
    """
    Perform a ncbi request.
    :param str id: Identifier of ncbi data set did you searching for.
    :param db: Related ncbi database to search into.
    :return: Returns request result.
    """
    result = []
    try:
        handle = Entrez.efetch(db=db, id=id, rettype='rb', retmode="xml")
        records = Entrez.parse(handle)
        for it in records:
            result.append(it)
        handle.close()
    except Exception as error:
        # TODO Logging system
        print(error, " on ", id, "\t", db)
    return result


def getByPath(obj, path):
    """
    Takes a ncbi response object and filter this by a given path.

    :param obj: Ncbi response object.
    :param str path: Path to a sub object oj the response object. Read about document explanation or documentation.
    :return: If path a valid path of obj the sub object will get returned. Otherwise, None.
    """
    if type(path) == str:
        path = path.split('.')
    if type(path) != list:
        return None

    for it in path:

        # Empty string mean take all
        if it == "":
            return obj

        # Dictionary iteration
        if it in obj:
            obj = obj[it]

        # List iteration
        elif it[0] == '[' and it[-1] == ']':
            try:
                if it.find('=') != -1:
                    temp = it[1:-1].split('=')[0:2]
                    for itt in obj:
                        if temp[0] in itt and itt[temp[0]] == temp[1]:
                            obj = itt
                            temp = False
                            break

                    # temp will set to False by success. If the loop finish without success the item could´t be found.
                    if temp:
                        return None
                else:
                    obj = obj[int(it[1:-1])]
            except:
                return None

        # Keyword search
        elif it[0] == '(' and it[-1] == ')':
            try:
                # Find one of this exclusive temp strings
                temp = it[1:-1].split(',')
                receive = None
                for itt in temp:
                    if json.dumps(obj).find(itt) != -1:

                        # We find a second name out of temp names.
                        if receive:
                            return None
                        else:
                            receive = itt
                return receive
            except:
                return None
        else:
            return None
    return obj


def printObj(obj, depth, spacer=""):
    """
    Command line prints of dictionaries are somtimes confusing, there for this function will print a dictionary
    like a tree list.

    :param obj: Any object like list and dictionaries.
    :param depth: Depth of the printed tree.
    :param spacer: Recursion param. Should never set by person.
    :return: void
    """
    if depth > 0:

        if type(obj) == dict or type(obj) == Entrez.Parser.DictionaryElement:
            for it in list(obj.keys()):
                print(spacer, it)
                printObj(obj[it], depth - 1, spacer + "\t")
        elif type(obj) == list or type(obj) == Entrez.Parser.ListElement:
            for it in obj:
                print(spacer, it)
                printObj(it, depth - 1, spacer + "\t")
        else:
            print(spacer, obj)
    else:
        return


if __name__ == '__main__':
    print("------ ncbi PATH Test -----")
    command = ""

    # Example notes ncbi559295 taxonomy, XP_033769018.1 protein

    while command.find("exit") == -1:
        command = input("\nChoose 'info','get','path', 'tree' ->")
        try:
            if command.find("info") != -1:
                print(getDbInfo(input("Database name ->")))

            elif command.find("get") != -1:
                res = get(input("Enter id ->"), input("Enter database ->"))
                print(res)

            elif command.find("path") != -1:
                res = get(input("Enter id ->"), input("Enter database ->"))
                print(res)
                while command.find("exit") == -1 or command.find("back") == -1:
                    print(getByPath(res[0], input("Enter path ->")))

            elif command.find("tree") != -1:
                res = get(input("Enter id ->"), input("Enter database ->"))
                print(res)
                printObj(res[0], int(input("Enter tree depht ->")))

            else:
                print("Invalid input. Close with 'exit'")


        except Exception as error:
            print(error)
