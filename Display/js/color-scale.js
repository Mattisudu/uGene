const COLOR_SCALE = {
  "sources": [
	{
		"Description":"INTRODUCTION TO COLOUR SCHEMES",
		"Author":  "Paul Tol (mail: p.j.j.tol@sron.nl)",
		"Link": "https://personal.sron.nl/~pault/#fig:scheme_bright"
	},
	{
		"Description":"Best Practices for Colour Blind Friendly Publications and Descriptions",
		"Author":  "Alba Fernández-Barral, CTAO Outreach and Education Coordinator",
		"Link": "https://www.cta-observatory.org/wp-content/uploads/2020/10/CTA_ColourBlindness_BestPractices-1.pdf"
	},
	{
		"Description":"Best series of colors to use for differentiating series in publication-quality plots",
		"Author":  "StackExchange.com Acount namens:  enigmaticPhysicist, Joseph D'Arimathea, a different ben",
		"Link": "https://stats.stackexchange.com/questions/118033/best-series-of-colors-to-use-for-differentiating-series-in-publication-quality"
	}
	
  ],
  "data": [
    {
      "label": "Rainbow",
      "value": "#DF0101 #FFFF00 #298A08 #00FF00 #01DFD7 #0101DF #F781BE"
    },
	{
      "label": "Mixed",
      "value": "#e6194b #3cb44b #ffe119 #0082c8 #f58231 #911eb4 #46f0f0 #f032e6 #d2f53c #fabebe #008080 #e6beff #aa6e28 #fffac8 #800000 #aaffc3 #808000 #ffd8b1 #000080 #808080 #ffffff #000000"
    },
    {
      "label": "Colorbrewer",
      "value": "#e41a1c #377eb8 #4daf4a #984ea3 #ff7f00 #ffff33 #a65628 #f781bf"
    },
	{
      "label": "Balance Eyesight",
      "value": "#77AADD #EE8866 #EEDD88 #FFAABB #99DDFF #44BB99 #BBCC33 #AAAA00 #DDDDDD"
    },
	{
      "label": "Bright Colorbilnd",
      "value": "#4477AA #EE6677 #228833 #CCBB44 #66CCEE #AA3377 #BBBBBB"
    },
	{
      "label": "Contrast Colorblind",
      "value": "#EE7733 #0077BB #33BBEE #EE3377 #CC3311 #009988 #BBBBBB"
    },
	{
      "label": "Vibrant Colorblind",
      "value": "#CC6677 #332288 #DDCC77 #117733 #88CCEE #882255 #44AA99 #999933 #AA4499"
    },
	{
      "label": "Silk Colorblind",
      "value": "#BBCCEE #CCEEFF #CCDDAA #EEEEBB #FFCCCC #DDDDDD #222255 #225555 #225522 #666633 #663333 #555555"
    },
	{
      "label": "Seaborn Colorblind",
      "value": "#0072b2 #009e73 #d55e00 #cc79a7 #f0e442 #56b4e9"
    },
	{
      "label": "Tableau Colorblind",
      "value": "#006ba4 #ff800f #ababab #595959 #5f9ed1 #c85200 #898989 #a2c8ec #ffbc79 #cfcfcf"
    },
	{
      "label": "Seaborn Dark",
      "value": "#1f77b4 #ff7f0e #2ca02c #d62728 #9467bd #8c564b #e377c2 #7f7f7f #bcbd23 #17becf"
    },
	{
      "label": "Bluesky",
      "value": "#2271b3 #89cff0 #3bbcd3 #058b8c"
    },
	{
      "label": "Colorcycle",
      "value": "#332288 #88CCEE #44AA99 #117733 #999933 #DDCC77 #CC6677 #882255 #AA4499"
    },
	{
      "label": "Orange Island",
      "value": "#E69F00 #56B4E9 #009E73 #F0E442 #0072B2 #D55E00 #CC79A7 #000000"
    }
  ]
}

function cVecAdd(c_vec_a, c_vec_b){
    /* Just add to arrays like vectors together.
    :return: Returns an array with addition values.*/
    try {
        return c_vec_a.map((it, i) => it +  c_vec_b[i])
    }catch (error){
        console.log(error)
        return Array.from(c_vec_a)
    }
}

function cVecSub(c_vec_a, c_vec_b){
    /* Just subtract to arrays like vectors each other.
    :return: Returns an array with c_vec_a - c_vec_b values.*/
    try {
        return c_vec_a.map((it, i) => it -  c_vec_b[i])
    }catch (error){
        console.log(error)
        return Array.from(c_vec_a)
    }
}

function cVecSkala(c_vec, scale){
    /* Just multiply an arrays with an scale.
    :return: Returns an array with c_vec_a * scale values.*/
    try {
        return c_vec.map((it, i) => it * scale)
    }catch (error){
        console.log(error)
        return Array.from(c_vec)
    }
}

// --------------  Clolor funcitons ----------------
function rgbStrToVec(color){
    /* Converts a hex color code string into a numpy 3 vector.
    :param color: Color code string. An example would be "#1A05FF".*/
    try{
    return [parseInt(color.slice(1,3), 16),
            parseInt(color.slice(3,5), 16),
            parseInt(color.slice(5,7), 16)]
    }catch (error){
        console.log(error)
        return [0, 0, 0]
    }
}

function rgbVecToStr(c_vec){
    /* Converts a numpy 3 vector within int variables into a rbg hex color code string.
    :return: Returns a string with a hex color code.*/
    try{

        if (c_vec[0] < 0)c_vec[0] = 0;
        if (c_vec[0] > 255)c_vec[0] = 255;

        if (c_vec[1] < 0)c_vec[1] = 0;
        if (c_vec[1] > 255)c_vec[1] = 255;

        if (c_vec[1] < 0)c_vec[1] = 0;
        if (c_vec[1] > 255)c_vec[1] = 255;

        return "#" + ("00"+ c_vec[0].toString(16)).slice(-2)+
                    ("00"+ c_vec[1].toString(16)).slice(-2) +
                    ("00"+ c_vec[2].toString(16)).slice(-2);
    }catch(error){
        console.log("Error: required_functionalities->rgbVecToStr()", error)
        return "#000000"
    }
}

function colorRampPalette(colors, n) {
    /*Interpolate colors linearly to create a color palette.
    :return:        Gives a list with hex color strings.*/
    result = []
    c_len = colors.length
    if (colors.length < 1) return [];
    if (colors.length == 1) return colors * n;
    if (n == 1) return [colors[0]];

    step = (colors.length - 1) / (n - 1)
    for (var i = 0; i < n; i++){
        if (Math.floor(step * i) == Math.ceil(step * i)) {
            result.push(colors[Math.floor(step * i)])
        }
        else{
            var v_color_a = rgbStrToVec(colors[Math.floor(step * i)]);
            var v_color_b = rgbStrToVec(colors[Math.ceil(step * i)])

            var v_color = cVecAdd(v_color_a , cVecSkala(cVecSub(v_color_b, v_color_a), (step * i % 1)))
            result.push(rgbVecToStr(v_color.map((it, i) => Math.round(it))))
        }
    }
    return result
}

function qualitativeColours(n, color_root=undefined){
    /* Generates a color palette in order to be able to differentiate between
    individual taxa as well as possible.
    :return:    Gives a list with hex color strings.*/
    const defauld_root = ["#DF0101", "#FFFF00", "#298A08", "#00FF00", "#01DFD7", "#0101DF", "#F781BE"]

    if (color_root === undefined){
        color_root = defauld_root
    }else {

        try {
            color_root = color_root.split()

            //Simply data check, not valid against none hex letters.
            for (i in color_root){
                if (color_root[i].length != 7 || color_root[i][0] != '#') {
                    console.log("Error: qualitiveColours() ", "ValueError: not a hex color string.")
                    color_root = defauld_root
                }}
        } catch (error) {
            //TODO Exchange print with common log function.
            console.log("Error: qualitiveColours() ", error)
            color_root = defauld_root
        }
    }
    return colorRampPalette(color_root, n)
}
