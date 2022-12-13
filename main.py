"""
This is the core file of uGene and performs the uMap clustering. All inputs being made by command line input.
Example: $> python outToCSV.py file=exampleData.out
"""

import pandas as pd
import umap
import json
import sys
import concurrent.futures

RESULT_COL_NAMES = []


def reduceList(li):
    """
    Unpack a list witch one contains list and join this list into one list.

    :param list li: Takes python standard list with in lists. [[a],[b]]
    :return: Returns a simply list. [a, b]"""
    result = []
    try:
        for it in li:
            result = result + it
    except Exception as error:
        # TODO Logging system
        print("Error reduceList()", error)
    return result


def jsonLoads(json_str, error_return):
    """
    Save use of json.load() to convert a json string to an object.

    :param str json_str: A simpy json string
    :param object error_return: A object witch one should get returned in a case of an error.
    :return: Give back the parsed json sting or the error_return object.
    """
    try:
        json_str = json_str.replace("'", '"')
        json_str = json.loads(json_str)
        return json_str
    except:
        # TODO Logging system
        print("Error can´t interpret ", json_str, " arguments.")
    # Error case gives origin string back.
    return error_return


def clusterDF(df, job_name='unnamedJob', angular_rp_forest=False, b=None,
              force_approximation_algorithm=False, init='spectral', learning_rate=1.0,
              local_connectivity=1.0, low_memory=False, metric='euclidean',
              metric_kwds=None, min_dist=0.1, n_components=2, n_epochs=None,
              n_neighbors=15, negative_sample_rate=5, output_metric='euclidean',
              output_metric_kwds=None, random_state=42, repulsion_strength=1.0,
              set_op_mix_ratio=1.0, spread=1.0, target_metric='categorical',
              target_metric_kwds=None, target_n_neighbors=-1, target_weight=0.5,
              transform_queue_size=4.0, transform_seed=42, unique=False, verbose=False):
    """
    Perform the main cluster job with the uMap cluster algorithm form umap-sklearn.

    :param df: Takes a pandas dataframe except only filled with raw data and organised like a matrix.
    :param str job_name: Just a name of this job. It is needed to create the result column name. Use one name just once.
    :param angular_rp_forest:
    :param b:
    :param force_approximation_algorithm:
    :param init:
    :param learning_rate:
    :param local_connectivity:
    :param low_memory:
    :param str metric: Means the metric witch one calculates the distance of data from each other. Ex.'euclidean'
    :param metric_kwds:
    :param float min_dist:
    :param int n_components: Defines how many dimensions this uMap should have.
    :param n_epochs:
    :param int n_neighbors: Heuristic value, witch effect a global vs local vision of the algorithm.
    :param negative_sample_rate:
    :param output_metric:
    :param output_metric_kwds:
    :param random_state:
    :param repulsion_strength:
    :param set_op_mix_ratio:
    :param spread:
    :param target_metric:
    :param target_metric_kwds:
    :param target_n_neighbors:
    :param target_weight:
    :param transform_queue_size:
    :param transform_seed:
    :param unique:
    :param verbose:
    :return:
    """
    filter_none_col = [it for it in RESULT_COL_NAMES if it in list(df.columns)]
    umap_result = umap.UMAP(
        angular_rp_forest=angular_rp_forest, b=b,
        force_approximation_algorithm=force_approximation_algorithm,
        init=init, learning_rate=learning_rate,
        local_connectivity=local_connectivity,
        low_memory=low_memory, metric=metric,
        metric_kwds=metric_kwds, min_dist=min_dist,
        n_components=n_components, n_epochs=n_epochs,
        n_neighbors=n_neighbors, negative_sample_rate=negative_sample_rate,
        output_metric=output_metric,
        output_metric_kwds=output_metric_kwds, random_state=random_state,
        repulsion_strength=repulsion_strength,
        set_op_mix_ratio=set_op_mix_ratio, spread=spread,
        target_metric=target_metric,
        target_metric_kwds=target_metric_kwds, target_n_neighbors=target_n_neighbors,
        target_weight=target_weight,
        transform_queue_size=transform_queue_size, transform_seed=transform_seed,
        unique=unique, verbose=verbose
    ).fit_transform(df.drop(filter_none_col, axis=1))

    # Save the UMAP output data
    cord = ['x', 'y', 'z']
    return [(job_name + str(n_components) + 'd_' + it[0], it[1]) for it in zip(cord, umap_result.transpose())]


def pivotDF(df, index, column, values, nan=float('nan'), report=None):
    """
    Manage the pandas.pivot() calculation. Allows column combinations and do conflict management. Index and column ar
    somthing like key in a dataframe. Conflicts appear whenever these keys are not unique.

    :param df: Standard molten dataframe.
    :param list index: List of column names witch one should affect the y-axis. Cluster subjects.
    :param list column: List of column names witch one should affect the x-axis.
    :param values: List with columns witch containing raw date for the matrix. All these columns should contain values.
    :param float nan: Value witch get filled into the matrix always there are no data available.
    :param bool report: Simple flac. If True, a simple report about conflict into the dataframe get printed.
    :return: Returns a pandas dataframe only filled with raw data and organised like a matrix.
    """
    #  Cut unimportant data columns
    df = df.drop([i for i in df.columns if i not in ([index, column] + values)], axis=1)

    # Report conflict deviation. Conflicts are given by the chosen key combination.
    # Would the key be a real key there are no conflicts will by definition.
    if report:
        groups_df = df.groupby([index, column])
        # Report max deviation, mean deviation and
        std_df = groups_df.std()
        size_series = groups_df.size()
        print("\n->Deviation report for database conflicts:")
        print("\tChosen keys : ", index, column)
        print("\tTotal count conflicts : ", str(size_series[size_series > 1].count()))
        print("\tMean of standard deviation : ", std_df.mean().to_string())
        print("\tTotal max deviation : ", std_df.max().to_string())
        df = groups_df.mean().reset_index()
    else:
        # Reduce Conflicts without do a conflict report.
        df = df.groupby([index, column]).mean().reset_index()

    if nan != float('nan'):
        return df.pivot(index=[index], columns=[column], values=values).fillna(nan)
    # else:
    return df.pivot(index=[index], columns=[column], values=values)


def mainAnalytics(df, x_axis=[], y_axis=[], values=[], jobs=[], dev_report=False, min_average=None):
    """
    This is a kind of main function witch one is called for every task. All over it effect program and processing
    structure.

    :param df: Standard molten pandas dataframe.
    :param str|list x_axis: One column name or list of column names witch one get used to pivot.
    :param str|list y_axis: One column name or list of column names witch one get used to pivot.
    :param str|list values: One column name or list of column names witch one get used to pivot.
    :param str|dict jobs: One sting witch effect the job_name or a job dictionary witch contains params of clusterDF()
    :param bool dev_report: Simple flac. If True, a simple report about conflict into the dataframe get printed.
    :param float min_average: Minimum average value of raw data into the data matrix to prevent a bad axis choice.
    :return: Returns a standard molten pandas dataframe with all results.
    """
    global RESULT_COL_NAMES
    # df,x_axis, y_axis, values -> How to interpret the data frame.
    # jobs -> [{job},{job_name: gen}]
    # job_name='tax', n_components=1, n_neighbors=15, min_dist=0.1, metric='euclidean'

    # Check inputs
    if len(x_axis) == 0 or len(y_axis) == 0 or len(values) == 0 or len(jobs) == 0:
        # TODO Logging system
        print("Error mainAnalytics() unable arguments x_axis, y_axis, values or jobs. ")
        return df

    # Process x_axis, y_axis and values given by a simple string. After these have to
    if type(x_axis) != list:
        x_axis = [str(x_axis)]
    if type(y_axis) != list:
        y_axis = [str(y_axis)]
    if type(values) != list:
        values = [str(values)]
    if type(jobs) != list:
        jobs = [{'job_name': str(jobs)}]

    # Columns are able to combine to get pivot keys.
    column = tuple(x_axis)
    index = tuple(y_axis)

    try:
        if len(x_axis) > 1:
            df[tuple(x_axis)] = df[x_axis].apply(lambda x: tuple(x), axis=1)
        else:
            column = x_axis[0]
            df[column]
        if len(y_axis) > 1:
            df[index] = df[y_axis].apply(lambda x: tuple(x), axis=1)
        else:
            index = y_axis[0]
            df[index]
    except Exception as error:
        # TODO Logging system
        print("Error invalid axis names. Axis have to be column names.\n", error)

    # The mat_df should only contain data values. All new cluster result column names will add to RESULT_COL_NAMES
    mat_df = pivotDF(df, column=column, index=index, values=values, nan=-1, report=dev_report)

    if dev_report or min_average is not None:
        sum_sum = mat_df.sum().sum()
        if dev_report:
            # If dev_report is ture, you may search for a good column pivot key, also you need a quality feedback.
            print("\n->Pivot matrix statistics:")
            print("\tCluster matrix size: x=", mat_df.shape[1], " y=", mat_df.shape[0])
            # Remember all NaN will set to -1, therefor a sum becomes a statement of the -1 fillment.
            print("\tCluster matrix sum: ", sum_sum)
            print("\tCluster matrix average value : ", sum_sum / mat_df.size)
            print("\tMatrix head:")
            print(mat_df.head())
        if min_average is not None and sum_sum / mat_df.size < float(min_average):
            # TODO logging system
            print("Abort task for axis ", index, " ", column, " by min_average border.")

            # CLean the dataframe from additional column and index. On single keys no additional column was created.
            if len(x_axis) > 1:
                df = df.drop([column], axis=1)
            if len(y_axis) > 1:
                df = df.drop([index], axis=1)

            return df

    # Do all the cluster jobs.

    res_df = pd.DataFrame(index=mat_df.index)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        res = [executor.submit(clusterDF, mat_df, **it) for it in jobs]

        # clusterDF returns a tuple with (column name, data list).
        for it in res:
            for itt in it.result():
                RESULT_COL_NAMES.append(itt[0])
                res_df[itt[0]] = itt[1]

    # Drop all column level. Old version bug fix and a case witch should never appear now.
    while res_df.columns.nlevels > 1:
        print("Error mainAnalytics()! Drop res_df multi level columns.")
        res_df = res_df.droplevel(1, axis=1)

    df = df.merge(res_df, left_on=[index], right_index=True, how="left")

    # CLean the dataframe from additional column and index. On single keys no additional column was created.
    # On change, make sure upper code duplicates will also change.
    if len(x_axis) > 1:
        df = df.drop([column], axis=1)
    if len(y_axis) > 1:
        df = df.drop([index], axis=1)

    return df


def main():
    # Default starting parameters.
    file_name = [""]
    tasks = ["[]"]
    # Argument pre processing
    pointer = [None]
    for i in range(1, len(sys.argv)):

        if sys.argv[i].find('file=') != -1:
            file_name[0] = sys.argv[i].replace('file=', '')
            pointer = file_name
        elif sys.argv[i].find('-f=') != -1:
            file_name[0] = sys.argv[i].replace('-f=', '')
            pointer = file_name
        elif sys.argv[i].find('mtask=') != -1:
            tasks[0] = sys.argv[i].replace('mtask=', '')
            pointer = tasks
        elif sys.argv[i].find('-t=') != -1:
            tasks[0] = sys.argv[i].replace('-t=', '')
            pointer = tasks
        else:
            pointer[0] = pointer[0] + sys.argv[i]

    # ----------------------- Debug ------------------------------------------
    # Bypass the command line input with the folowing lines.
    # job_name='unnamedJob', n_components=1, n_neighbors=15, min_dist=0.1, metric='euclidean'
    #file_name = "test.csv"
    #tasks = [{'dev_report': 1, 'x_axis': 'ncbiID', 'y_axis': 'geneID', 'values': 'FAS_F',
    #          'jobs': [{'job_name': 'gene'}, {'job_name': 'gene', 'n_components': 1},{'job_name': 'gene', 'n_components': 3},
    #                   {'job_name': 'geneCor', 'metric':'correlation'}, {'job_name': 'geneCor', 'n_components': 1, 'metric':'correlation'},{'job_name': 'geneCor', 'n_components': 3, 'metric':'correlation'}]
    #          }]
    #

    # Argument post processing
    file_name = str(file_name[0])
    tasks = jsonLoads(tasks[0], [])
    if not tasks:
        return

    # Analytics processing
    try:
        print("Load Data ", file_name)
        df = pd.read_csv(file_name)

        # Print all dataframe columns.
        # pd.set_option('display.max_columns', None)
        print(df.head())
    except Exception as error:
        # TODO Logging system
        print("Cant´s open '", file_name, "' file.\n", error)
        return;

    for cur_task in tasks:
        try:
            print("Process task : ", cur_task)
            df = mainAnalytics(df, **cur_task)

        except Exception as error:
            print("Task fail! \t", error)

    # print(df.columns)
    # print(df.head)

    df.to_csv(file_name.replace(".csv", ".cluster.csv"))


if __name__ == "__main__":
    main()
