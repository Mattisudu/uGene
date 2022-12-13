var DataFrame = dfjs.DataFrame;
function main (){

    //This tool is´t ready done jet.
    // To jused anyway, edit file name and take care that the loaded datefaem contains the used column names.
    initOutputTable('select_output')
    mainInitPlott('plott_1_container', 'exapleData.cluster.csv','geneID','geneID', 'gene2d_x', 'gene2d_y', 'select_output')
    mainInitPlott('plott_2_container', 'exapleData.cluster.csv','geneID','geneID', 'gene2d_x', 'gene2d_y', 'select_output')
}

function progressbarSetOrg(id = "", progress = undefined){
    const elem = document.getElementById(id);
    console.log(elem)
    console.log(elem.querySelector("#sub-first"))
    if (elem !== null && elem.querySelector("#sub-first") !== null) {

        if (progress === undefined || progress === NaN || progress < 0) {
            elem.classList.add("hidden");
            elem.classList.remove('progress')
            elem.querySelector("#sub-first").innerText = ""
        } else {
            elem.classList.remove("hidden")
            elem.classList.add('progress')
            elem.querySelector("#sub-first").style.width = String(progress) + "%"
            elem.querySelector("#sub-first").innerText = String(progress) + "%"
        }
    }
}

function progressbarSet(elem, progress = undefined){

    if (elem !== null && elem.querySelector("._progress_bar_") !== null) {

        if (progress === undefined || progress === NaN || progress < 0) {
            elem.classList.add("hidden");
            elem.classList.remove('progress')
            elem.querySelector("._progress_bar_").style.width = "0%"
            elem.querySelector("._progress_bar_").innerText = ""
        } else {
            elem.classList.remove("hidden")
            elem.classList.add('progress')
            elem.querySelector("._progress_bar_").style.width = String(progress) + "%"
            elem.querySelector("._progress_bar_").innerText = String(progress) + "%"
        }
    }
}

function mainInitPlott(plott_container,csv_file, col_main, col_color, col_data_x, col_data_y, output_table = undefined){
    const container = document.getElementById(plott_container)
    if (container === null || container === undefined) return false;

    const plott_id = container.querySelector('._plott_').getAttribute('id')
    const pro = container.querySelector('._progress_')

    if (plott_id === null || pro === null) return false;

    console.log("PlotID:", plott_id)

    const init = initPlot(plott_id, csv_file, col_main, col_color, col_data_x, col_data_y,pro, output_table)

    // After the plot is created we want to add some functions.
    init.then(()=>{
        initPlottDropdown(plott_container)
    });
}

function initPlottDropdown(plott_container){
    const container = document.getElementById(plott_container)
    const plot = container.querySelector('._plott_');
    const plot_id = plot.getAttribute('id')

    if (container === null || container === undefined) return false;

    if (plot.getAttribute('p_df') !== undefined && DATAFRAMES[plot.getAttribute('p_df')].listColumns() !== undefined) {
        const arr_col = Array.from(DATAFRAMES[plot.getAttribute('p_df')].listColumns())
        const new_dd = arr_col.map((it) =>{return "<li><button type=\"button\" class=\"dropdown-item\" onClick=\"updateColor('" + String(plot_id) + "','" + String(it) + "')\">" + String(it) + "</button></li>"}).join('');
        container.querySelector("._dd_color_ul_").innerHTML = new_dd
    }
}

function setOutputTable(table_id, arr_head, arr_rows){
    const table = document.getElementById(table_id)

    const new_head = "<tr>" + arr_head.map((it)=> "<th scope='col'>"+ String(it) +"</th>").join('') + "</tr>"
    const new_body = arr_rows.reduce((row_a, row_b)=>{
        return row_a + "<tr>" + row_b.map((it)=> "<th scope='col'>"+ String(it) +"</th>").join('') + "</tr>"
    },"")

    table.querySelector('._thead_').innerHTML= new_head
    table.querySelector('._tbody_').innerHTML= new_body
}

function initOutputTable(table){
    // Try to fix an id given plot.
    if(typeof table == "string"){
        table = document.getElementById(table)
    }

    // Checks if the plot div have the data.
    if (table === undefined) {
        console.log("initOutputTable() fail by id!")
        return false;
    }

    if (! table.getAttributeNames().includes('col_main')){
        const attr = document.createAttribute('col_main')
        attr.value = undefined
        table.setAttributeNode(attr)
    }
    if (! table.getAttributeNames().includes('list_main')){
        const attr = document.createAttribute('list_main')
        attr.value = undefined
        table.setAttributeNode(attr)
    }
}
function downloadOutputTable(table, filename = "table.csv"){// Try to fix an id given plot.
    if(typeof table == "string"){
        table = document.getElementById(table)
    }

    // Checks if the plot div have the data.
    if (table === undefined) {
        console.log("initOutputTable() fail by id!")
        return false;
    }

    const col_main = table.getAttribute('col_main')
    const main_list = table.getAttribute('list_main').split(',')
    const df = DATAFRAMES[0]

    if(df === undefined || main_list === undefined|| main_list === [] || col_main === undefined) return false;
    var down = document.createElement('a');
    down.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(df.filter(row => main_list.includes(row.get(col_main))).toCSV(true)));
    down.setAttribute('download', filename);

    down.style.display = 'none';
    document.body.appendChild(down);
    down.click();
    document.body.removeChild(down);

    return true
}
