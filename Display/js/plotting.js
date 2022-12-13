var DataFrame = dfjs.DataFrame;
var DATAFRAMES = []

function  plotting(){
    var xArray = [50,60,70,80,90,100,110,120,130,140,150];
    var yArray = [7,8,8,9,9,9,10,11,14,14,15];

// Define Data
    var data = [{
        x: xArray,
        y: yArray,
        mode:"markers",
        type:"scatter"
    }];

// Define Layout
    var layout = {
        xaxis: {range: [40, 160], title: "Square Meters"},
        yaxis: {range: [5, 16], title: "Price in Millions"},
        title: "House Prices vs. Size"
    };

    Plotly.newPlot("myPlot", data, layout);
}

function debugPlott(){
    return DataFrame.fromCSV('mini_test.cluster.csv').then((df) =>{

        df  = reduceDF2(df, 'geneID')
        //  -----  Load job is done. This reduce is for the main cluster axis.  ------

        const color_col = 'superkingdom'
        const traces = df.groupBy(color_col).toCollection()

        const color = qualitativeColours(traces.length)
        const data = traces.map((trace , i) => newTrace(trace,'gene2d_x', 'gene2d_y', color[i]))
        //const data = df.toCollection().map((t, i) => newTrace(t,'gene2d_x', 'gene2d_y', 'superkingdom', color[i]))
        const layout = newLayout('gene2d_x', 'gene2d_y')

        console.log("Plott:")
        console.log(data)
        console.log(layout)
        Plotly.newPlot("debugPlott", data, layout);
    })}


var helpReduce = function (arr_a, arr_b){arr_a.map((it, i) => it.add(arr_b[i])); return arr_a}
function updateHelpReduce(n_arr_length){
    let add = "";
    for(let i = 0; i < n_arr_length; i++){
        add = add + "arr_a["+String(i)+"].add(arr_b["+String(i)+"]);"
    }
    helpReduce = eval('(function(arr_a, arr_b) {'+add+'return arr_a;})')
}
function reduceDF(df, key = 'gevneID' , max_col = 30, pro = undefined) {
    if (pro !== undefined && pro !== null) progressbarSet(pro,0);

    // Map all columns at once to reduce runtime issues by map()

    const groups = df.groupBy(key).toCollection()

    let pro_value = 0;
    let pro_step = 100 / groups.length;

    var mapped = groups.map(function (_ref3) {
        const keys = _ref3.group.listColumns()
        const arr_length  = keys.length
        updateHelpReduce(arr_length)
        const init = Array.from(Array(arr_length).keys()).map((it)=> new Set())
        let res = _ref3.group.toArray().reduce(helpReduce,init)

        if(pro === undefined || pro === null){
            res = res.map((item_s, i) => {
                let result = ""
                for (const item of item_s){
                    result+=", " + String(item);
                    if(result.length >= max_col){
                        result+= " ... ("+String(item_s.size)+")";
                        break;
                    }
                }
                return result.substring(2);
            })
        }else{

            res = res.map((item_s, i) => {
                let result = ""
                for (const item of item_s){
                    result+=", " + String(item);
                    if(result.length >= max_col){
                        result+= " ... ("+String(item_s.size)+")";
                        break;
                    }
                }
                return result.substring(2);
            })
            progressbarSet(pro, pro_value)
            pro_value = pro_value + pro_step;
        }
        return res;
    });

    if (pro !== undefined && pro !== null) progressbarSet(pro, 100);
    return df.__newInstance__(mapped, df.listColumns());
}
function reduceDfSort(df, key = 'gevneID' , max_col = 30, pro = undefined) {
    if (pro !== undefined && pro !== null) progressbarSet(pro,0);

    // Map all columns at once to reduce runtime issues by map()

    const groups = df.groupBy(key).toCollection()

    let pro_value = 0;
    let pro_step = 100 / groups.length;

    var mapped = groups.map(function (_ref3) {
        const keys = _ref3.group.listColumns()
        const arr_length  = keys.length
        updateHelpReduce(arr_length)
        const init = Array.from(Array(arr_length).keys()).map((it)=> new Set())
        let res = _ref3.group.toArray().reduce(helpReduce,init)

        if(pro === undefined || pro === null){
            res = res.map((item_s, i) => {
                let result = ""
                item_s = Array.from(item_s).sort()
                for (const item of item_s){
                    result+=", " + String(item);
                    if(result.length >= max_col){
                        result+= " ... ("+String(item_s.length)+")";
                        break;
                    }
                }
                return result.substring(2);
            })
        }else{

            res = res.map((item_s, i) => {
                let result = ""
                item_s = Array.from(item_s).sort()
                for (const item of item_s){
                    result+=", " + String(item);
                    if(result.length >= max_col){
                        result+= " ... ("+String(item_s.length)+")";
                        break;
                    }
                }
                return result.substring(2);
            })
            progressbarSet(pro, pro_value)
            pro_value = pro_value + pro_step;
        }
        return res;
    });

    if (pro !== undefined && pro !== null) progressbarSet(pro, 100);
    return df.__newInstance__(mapped, df.listColumns());
}
function newTrace(trace, x, y, marker_color){

    const data = trace.group.toDict()
    const hovertemplate = Object.keys(data).map((key, i) => key +": %{customdata["+String(i) +"]}").join("<br>") + "<extra></extra>"

    //console.log("Trace:")
    //console.log(data)
    return {
                'customdata': trace.group.toArray(),
                'hovertemplate': hovertemplate,
                'legendgroup': Object.values(trace.groupKey)[0],
                'marker': {'color': marker_color, 'symbol': 'circle'},
                'mode': 'markers',
                'name': Object.values(trace.groupKey)[0],
                'orientation': 'v',
                'showlegend': true,
                'type': 'scatter',
                'x': data[x],
                'xaxis': 'x',
                'y': data[y],
                'yaxis': 'y'
            }
}
function newLayout(name_xaxis="", name_yaxis = ""){
    return {'legend': {'title': {'text': 'Legend'}, 'tracegroupgap': 0},
            'template': '...',
            'title': {'text': 'unnamedJob  '},
            'xaxis': {'anchor': 'y', 'domain': [0.0, 1.0], 'title': {'text': name_xaxis}},
            'yaxis': {'anchor': 'x', 'domain': [0.0, 1.0], 'title': {'text': name_yaxis}}}
}

function initPlot(plot_id,csv_file, col_main, col_color, col_data_x, col_data_y, pro = undefined, output_table = undefined){

    const plot = document.getElementById(plot_id)

    // Checks if the plot div have the data.
    if (plot === undefined) {
        console.log("initPlot() fail by plot_id!")
        return false;
    }

    if (!initPlotData(plot)) return false;
    // Attention if you use different filenames.
    return DataFrame.fromCSV(csv_file).then((df) =>{

        // Use a array with Dataframes to reduce reads
        if(DATAFRAMES[0] === undefined){
            DATAFRAMES[0] = df
        }
        df  = reduceDfSort(df, col_main, 30, pro)
        //  -----  Load job is done. This reduce is for the main cluster axis.  ------

        const p_DATAFRAMES = DATAFRAMES.length
        DATAFRAMES.push(df)

        plot.setAttribute('p_df', p_DATAFRAMES)

        plot.setAttribute('col_main', col_main)
        plot.setAttribute('col_color', col_color)

        plot.setAttribute('col_data_x', col_data_x)
        plot.setAttribute('col_data_y', col_data_y)

        console.log("set output_table", output_table)
        plot.setAttribute('output_table', output_table)

        plot.parentNode.parentNode.querySelector("._dd_main_").innerHTML = "Dataframe (" + col_main + ")"
        plot.parentNode.parentNode.querySelector("._dd_color_").innerHTML = "Color (" + col_color + ")"


        const traces = df.groupBy(col_color).toCollection()

        const color = qualitativeColours(traces.length)

        const data = traces.map((trace , i) => newTrace(trace,col_data_x, col_data_y, color[i]))
        //const data = df.toCollection().map((t, i) => newTrace(t,'gene2d_x', 'gene2d_y', 'superkingdom', color[i]))
        const layout = newLayout(col_data_x, col_data_y)

        console.log("Plott:")
        console.log(data)
        console.log(layout)
        Plotly.newPlot(plot_id, data, layout);
        registerPlotEvents(plot_id)

        if(pro !== undefined && pro !== null) progressbarSet(pro, -1);
        return df.listColumns();
    })

}
function updateMain(plot_id = "", new_main  = "", col_data_x = "", col_data_y = "", max_col = 30, pro=undefined){
    const plot = document.getElementById(plot_id);
    if(plot === undefined || plot === null) return;
    if (!checkPlotData(plot) || DATAFRAMES[0] === undefined) return;

    if(pro === undefined || pro === null){
        pro = plot.parentNode.parentNode.querySelector('._progress_')
    }

    plot.setAttribute('col_data_x', col_data_x)
    plot.setAttribute('col_data_Y', col_data_y)
    const col_color = plot.getAttribute('col_color')

    // Only set new col_main if it is required, to reduce calculation consumption.
    if (plot.getAttribute('col_main') !== new_main) {
        plot.setAttribute('col_main', new_main)

        // Delete all dataframe form global DATAFRAME storage.
        const old_p_DATAFRAME = plot.getAttribute('p_df')
        DATAFRAMES.splice(old_p_DATAFRAME, 1)

        // Build dataframe new with changed main base.
        DATAFRAMES.push(reduceDfSort(DATAFRAMES[0], new_main, max_col, pro))
        plot.setAttribute('p_df', String(DATAFRAMES.length - 1))

        // Set new text into dropdown main open button.
        plot.parentNode.parentNode.querySelector("._dd_main_").innerHTML = "Dataframe (" + new_main + ")"
    }

    // Set the old color again to refresh the plot.
    updateColor(plot_id,col_color)
    if(pro !== undefined && pro !== null) progressbarSet(pro, -1);
}

function updateColor(plot_id = "", new_color  = ""){
    const plot = document.getElementById(plot_id);
    if (!checkPlotData(plot)) return;


    plot.setAttribute('col_color', new_color)

    const col_data_x = plot.getAttribute('col_data_x')
    const col_data_y = plot.getAttribute('col_data_Y')
    let df = DATAFRAMES[parseInt(plot.getAttribute('p_df'))]

    // Set new text into dropdown open button
    plot.parentNode.parentNode.querySelector("._dd_color_").innerHTML = "Color (" + new_color + ")"

    const traces = df.groupBy(new_color).toCollection()
    const color = qualitativeColours(traces.length)

    const data = traces.map((trace , i) => newTrace(trace,col_data_x, col_data_y, color[i]))
    const layout = newLayout(col_data_x, col_data_y)

    console.log("Updata Color()")
    Plotly.newPlot(plot_id, data, layout);
    registerPlotEvents(plot)

    return df.listColumns()
}
function initPlotData(plot){
    // Try to fix an id given plot.
    if(typeof plot == "string"){
        plot = document.getElementById(plot)
    }

    // Checks if the plot div have the data.
    if (plot === undefined) {
        console.log("AddPlotData() fail by id!")
        return false;
    }

    if (! plot.getAttributeNames().includes('p_df')){
        const attr = document.createAttribute('p_df')
        attr.value = undefined
        plot.setAttributeNode(attr)
    }
    if (! plot.getAttributeNames().includes('col_main')){
        const attr = document.createAttribute('col_main')
        attr.value = undefined
        plot.setAttributeNode(attr)
    }
    if (! plot.getAttributeNames().includes('col_color')){
        const attr = document.createAttribute('col_color')
        attr.value = undefined
        plot.setAttributeNode(attr)
    }
    if (! plot.getAttributeNames().includes('col_data_x')){
        const attr = document.createAttribute('col_data_x')
        attr.value = undefined
        plot.setAttributeNode(attr)
    }
    if (! plot.getAttributeNames().includes('col_data_y')){
        const attr = document.createAttribute('col_data_y')
        attr.value = undefined
        plot.setAttributeNode(attr)
    }
    if (! plot.getAttributeNames().includes('output_table')){
        const attr = document.createAttribute('output_table')
        attr.value = undefined
        plot.setAttributeNode(attr)
    }

    return true
}
function checkPlotData(plot){

    // Try to fix an id given plot.
    if(typeof plot == "string"){
        plot = document.getElementById(plot)
    }

    // Checks if the plot div have the data.
    if (plot === undefined) {
        console.log("Fail to find plot by id !")
        return false;
    }

    // Check all Attributes by this own
    if ( !plot.getAttributeNames().includes('p_df') || !plot.getAttributeNames().includes('col_main') || !plot.getAttributeNames().includes('col_color')|| !plot.getAttributeNames().includes('col_data_x')|| !plot.getAttributeNames().includes('col_data_y')|| !plot.getAttributeNames().includes('output_table')){
        addPlotData(plot)
        return false
    }

    // output_table can be undefined....
    if( plot.getAttribute('p_df') === undefined || plot.getAttribute('col_main') === undefined || plot.getAttribute('col_color') === undefined || plot.getAttribute('col_data_x') === undefined || plot.getAttribute('col_data_y') === undefined){
        return false
    }
    return true
}

function registerPlotEvents(plot){
    // Try to fix an id given plot.
    if(typeof plot == "string"){
        plot = document.getElementById(plot)
    }

    // Checks if the plot div have the data.
    if (plot === undefined) {
        console.log("Fail to find plot by id !")
        return false;
    }
    const p_main = DATAFRAMES[plot.getAttribute('p_df')].listColumns().indexOf(plot.getAttribute('col_main')) | 0
    document.getElementById(plot.getAttribute('output_table')).setAttribute('col_main', plot.getAttribute('col_main'))


    plot.on('plotly_click', function(data){
        if(plot.getAttribute('output_table') !== undefined) {
            const df = DATAFRAMES[plot.getAttribute('p_df')]
            const col = df.listColumns()

            let data_list = []
            let main_list = []
            for (var i = 0; i < data.points.length; i++) {

                data_list.push(data.points[i].customdata)
                main_list.push(data.points[i].customdata[p_main])
            }

            // Show clicked data on table below
            setOutputTable(plot.getAttribute('output_table'), col, data_list)
            document.getElementById(plot.getAttribute('output_table')).setAttribute('list_main', main_list)
        }
    });
    plot.on('plotly_selected', function(data){
        if(plot.getAttribute('output_table') !== undefined) {
            const df = DATAFRAMES[plot.getAttribute('p_df')]
            const col = df.listColumns()

            let data_list = []
            let main_list = []
            for (var i = 0; i < data.points.length; i++) {
                data_list.push(data.points[i].customdata)
                main_list.push(data.points[i].customdata[p_main])
            }

            // Show clicked data on table below
            setOutputTable(plot.getAttribute('output_table'), col, data_list)
            document.getElementById(plot.getAttribute('output_table')).setAttribute('list_main', main_list)
        }
    });

    return true;
}

