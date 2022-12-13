"""
This is the fist part of the uGene pipline and is used to convert phylorbf .out files to standard csv files.
Additionally, this script is able to collect all kind of date from ncbi or merge additional data from other .csv files.

Example: $> python outToCSV.py file=exampleData.out -ncbi=[{'col':'ncbiID','db':'taxonomy','path':'lineage'}]
Example: $> python outToCSV.py file=test.out -cd=\u007C -e={'orthoID':['','temp','protID']}
"""
import sys
import pandas as pd
import json
import os
from Bio import Entrez


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


def ncbiLoad(df, col, db, path, mail="max@muellmail.com", rename=None):
    """
    Downloads data from ncbi.
    :param df: Pandas molten dataframe were the downloaded date should be added.
    :param col: Column of df witch one contains id´s witch should be searched for on ncbi.
    :param db: Related ncbi database to search into.
    :param str|list path: One or serval filters witch data from the ncbi response object being kept. Get by ncbiTest.py.
    :param mail: Your mail address if you want to get into contact if something went wrong.
    :param str|list rename: All new added columns get named by path, in order of paths rename them with this param.
    :return: Returns the df with all downloaded new filtered data.
    """
    try:
        if type(path) == str:
            path = [path]
        if type(rename) == str:
            rename = [rename]

        Entrez.email = mail
        id_list = df[col].unique().tolist()

        handle = Entrez.efetch(db=db, id=",".join(id_list), rettype='rb', remode='xml')

        records = Entrez.parse(handle)

        result = {i: [] for i in path}
        for it in records:
            for itt in path:
                # it is the incoming data. itt is path or filter of this data.
                result[itt].append(getByPath(it, itt))
        handle.close()

        if len(result[path[0]]) != len(id_list):
            # TODO find out how it could happend and how to prevent. This case should never appear.
            print("Bad request, no consistent ncbi response!")
            print("len(response) != len(requests) \t(", len(result[path[0]]), ") != (", len(id_list), ")")
            return df

        # In some cases we will not have the path equals the column name, there for we set to rename string.
        if rename:
            for it in list(zip(path, rename)):
                result[it[1]] = result.pop(it[0])

        # Add the merge id´s to the result dict.
        result[col] = id_list

        result = pd.DataFrame(result)

        return df.merge(result, left_on=col, right_on=col)
    except Exception as error:
        # TODO Logging system
        print("Fail ncbiLoad() ", error)

        # Bug fix, wrong error report
        if str(error).find("file should be opened in binary mode") != -1:
            print("This error accrue if the Ncbi Database is´t correct. Check this out!")

        return df


def importCSV(df, file="", col="", imports=[]):
    """
    Merge data out of a different csv file into the current working dataframe.

    :param df: Standard pandas molten dataframe.
    :param str file: File name of the csv file with the additional data.
    :param str col: Column name of these column witch contains the keys to merge the data.
    :param str|list imports: Just one or a list of column names witch get imported.
    :return: Returns the dataframe with the new additional data.
    """
    if type(imports) == str:
        imports = [imports]

    df_import = pd.read_csv(file)
    if not col in df_import.columns:
        print("Fail merge dataframes!")
        return df

    if imports == []:
        pass
    elif col in imports:
        df_import = df_import[imports];
    else:
        df_import = df_import[[col] + imports]

    return df.merge(df_import, left_on=col, right_on=col, how='left')


def main():
    """
    :arg edit: Gives a possibility to split up columns into several columns.
    :arg command: Contains the working file name and other program leading stings.
    :arg delimiter: Is just a char witch is the old delimiter char.
    :arg col_delimiter: Is just a char witch is the old delimiter in edit able columns.
    :return: void
    """
    ncbi_pre_load = []
    ncbi_load = []
    edit = {}
    command = "exit"
    delimiter = '\t'
    col_delimiter = '|'
    merge = []

    # Argument pre processing
    for i in range(1, len(sys.argv)):
        print(sys.argv[i], end=" ")

        if sys.argv[i].find('file=') != -1:
            command = sys.argv[i].replace('file=', '')
        if sys.argv[i].find('-f=') != -1:
            command = sys.argv[i].replace('-f=', '')

        if sys.argv[i].find('edit=') != -1:
            edit = sys.argv[i].replace('edit=', '')
        if sys.argv[i].find('-e=') != -1:
            edit = sys.argv[i].replace('-e=', '')

        if sys.argv[i].find('-d=') != -1:
            delimiter = sys.argv[i].replace('-d=', '')
        if sys.argv[i].find('delimiter=') != -1:
            delimiter = sys.argv[i].replace('-delimiter=', '')

        if sys.argv[i].find('-cd=') != -1:
            col_delimiter = sys.argv[i].replace('-cd=', '')
        if sys.argv[i].find('cdelimiter=') != -1:
            col_delimiter = sys.argv[i].replace('cdelimiter=', '')

        if sys.argv[i].find('-ncbi=') != -1:
            ncbi_load = sys.argv[i].replace('-ncbi=', '')
        if sys.argv[i].find('ncbiload=') != -1:
            ncbi_load = sys.argv[i].replace('ncbiload=', '')

        if sys.argv[i].find('-prencbi=') != -1:
            ncbi_pre_load = sys.argv[i].replace('-prencbi=', '')
        if sys.argv[i].find('prencbiload=') != -1:
            ncbi_pre_load = sys.argv[i].replace('prencbiload=', '')

        if sys.argv[i].find('-m=') != -1:
            merge = sys.argv[i].replace('-m=', '')
        if sys.argv[i].find('merge=') != -1:
            merge = sys.argv[i].replace('merge=', '')

    # Argument post processing
    if edit != {}:
        try:
            edit = edit.replace("'", '"')
            edit = json.loads(edit)
        except Exception as error:
            edit = {}
            # TODO Logging system
            print("Error can´t read edit arguments.\n", error)

    if ncbi_load != []:
        try:
            ncbi_load = json.loads(ncbi_load.replace("'", '"'))
        except Exception as error:
            ncbi_load = []
            # TODO Logging system
            print("Error can´t read ncbi arguments.\n", error)

    if ncbi_pre_load != []:
        try:
            ncbi_pre_load = json.loads(ncbi_pre_load.replace("'", '"'))
        except Exception as error:
            ncbi_pre_load = []
            # TODO Logging system
            print("Error can´t read ncbi arguments.\n", error)

    if merge != []:
        try:
            merge = json.loads(merge.replace("'", '"'))
        except Exception as error:
            merge = []
            # TODO Logging system
            print("Error can´t read ncbi arguments.\n", error)

    # Delimiter have to be just one char
    delimiter = delimiter[0]
    # TODO Fix for special chars they can´t given in comand line. Feel free to add more unicode translation.
    col_delimiter = str(col_delimiter).replace('\\u007C', '|')
    col_delimiter = str(col_delimiter).replace('u007C', '|')
    col_delimiter = list(col_delimiter)

    print("Command: ", command)
    print("Edit: ", edit)
    print("Delimiter: ", delimiter)
    print("ColumnDelemitter: ", col_delimiter)
    print("Pre NCBI load: ", ncbi_pre_load)
    print("NCBI load: ", ncbi_load)
    print("Merge:", merge)

    if command == "exit":
        # TODO Edit log system
        print('Exit by no command')
        return;

    # Translation part
    try:
        if command.find(".out") != -1:
            with open(command, 'r') as file:
                data = file.read().replace(delimiter, ',')
                # TODO Add a none .out  to csv converting mode.

                # Ask for edit calculation ? Without no .temp file will be created.
                if len(edit) > 0:
                    command = command.replace(".out", ".temp.csv")
                else:
                    command = command.replace(".out", ".csv")
                with open(command, 'w') as file2:
                    file2.write(data)

        # Edit calculation part.
        if command.find(".csv") != -1 and (
                len(edit) > 0 or len(ncbi_load) > 0 or len(ncbi_pre_load) > 0 or len(merge) > 0):
            df = pd.read_csv(command)

            if "Unnamed: 0" in df.columns:
                df = df.drop("Unnamed: 0", axis=1)

            # Do the preload ncbi jobs
            for it in ncbi_pre_load:
                df = ncbiLoad(df, **it)

            # Index to current col_delimiter
            p_col_delimiter = 0
            for it in edit:
                print("Edit column: ", it)

                if p_col_delimiter >= len(col_delimiter):
                    p_col_delimiter = len(col_delimiter) - 1
                # it contains the line to need to edit
                df[it] = df[it].str.split(col_delimiter[p_col_delimiter])

                # Increment now. Don´t use p_col_delimiter for the rest of the loop.
                p_col_delimiter = p_col_delimiter + 1

                for itt in range(0, len(edit[it])):
                    # Space or empty column names are not aloud.
                    if edit[it][itt] != "" and edit[it][itt] != "" and edit[it][itt] != it:
                        # itt contains the iterated list index of edit[it]
                        df[edit[it][itt]] = df[it].apply(lambda x: x[itt])

                # Process column name override.
                if it in edit[it] and it.strip():
                    print("Column override")
                    itt = edit[it].index(it)
                    df[edit[it][itt]] = df[it].apply(lambda x: x[itt])
                else:

                    # Drop redundant columns.
                    df = df.drop(it, axis=1)

            # Do some merge stuff
            for it in merge:
                df = importCSV(df, **it)

            # Do the normal ncbi jobs
            for it in ncbi_load:
                df = ncbiLoad(df, **it)

            if command.find('.temp.csv') != -1:
                # Remove .temp file.
                try:
                    os.remove(command)
                except Exception as error:
                    # TODO Logging system
                    print(error)

                df.to_csv(command.replace(".temp", ""), index=False)
            else:
                # There is no temp file because may there was no .out file convention.
                df.to_csv(command, index=False)

    except Exception as e:
        # TODO Edit Logging system
        print(e)

    print("outToCSV.py ... Finish!")


main()
