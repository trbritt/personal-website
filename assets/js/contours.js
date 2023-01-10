const queryString = window.location.search;
const urlParams = new URLSearchParams(queryString);
var stringSeed = "tristan";
if(urlParams.has('seed')) {
    stringSeed = urlParams.get('seed');
    console.log('Using seed: ' + stringSeed);
}
let pixeld = 1;
var stringSeedHash = cyrb128(stringSeed);
var fxrand = sfc32(...stringSeedHash);
let callpreviewtimer;
let colorsarray  = [
    [["1"] , ["2b2d42","8d99ae","ef233c","d90429"]],
    [["2"] , ["0A3A4A","196674","33A6B2","9AC836"]],
    [["3"] , ["91D9CC","F2C48D","593520"]],
    [["4"] , ["F24968","D984A3","BDDEF2","DCEAF2","A61212"]],
    [["5"] , ["BAD9D9","7BA692","344023","AFBF34"]],
    [["6"] , ["A62940","D90B42","591E3A","D8D9C5","D98484"]],
    [["7"] , ["F2C438","F27B13","D9C7B8","8C6C5A","F26D3D"]],
    [["8"] , ["F20544","8C2656","3F618C","012340","011826"]],
    [["9"] , ["BF567D","8C4D70","435C73","6FA8BF","539DA6"]],
    [["10"] , ["F2711B","B7C093","97A275"]],
    [["11"] , ["593E40","F2D5C4","F2785C","F2594B","A66F6F"]],
    [["12"] , ["D8E3CF","F0A36F","C5557B","4F2C62"]],
    [["13"] , ["1C2747","EC6508","1375BB","E9483B","D8DDE5"]],
    [["14"] , ["8C2354","F2CB05","F29F05","F27405"]],
    [["15"] , ["64C5F5","2484BF","4DA60D","2A5908"]],
    [["16"] , ["04ADBF","04BFBF","025959","A0A603"]],
    [["17"] , ["5B98A6","026873","32D9D9","0D0D0D"]],
    [["18"] , ["26261C","595843","A6A485","0D0D0D"]],
    [["19"] , ["3D5A73","4F6F8C","F2E2CE","F24130","8C0303"]],
    [["20"] , ["111440","144E73","1F80A6","469CA6","F2F2F2"]],
    [["21"] , ["26261C","595843","A6A485","0D0D0D"]],
    [["22"] , ["4C291E","FFB632","D2430C","872018"]],
    [["23"] , ["8AB0BF","61888C","17261E","3B593F","467339"]],
    [["24"] , ["592040","733462","2A708C","50ABBF","0D0D0D"]],
]
let colorarray = randomarray(colorsarray);
let colors = colorarray[1];
let colorsshuffle = shufflearr(colors);
let pointarray = [];
let bgc ;
let _noiseCounter = 9999;
let _noiseFloor = 0.1;
let _width = 2000;
let _height = 2800;
let _singleWidth = _singleHeight = rndint(15, 60);
let _singleWidthDiv = 1/_singleWidth;
let _minDist = rndint(_singleWidth*0.4, _singleWidth*1.8);
let _maxLen = rndint(_singleWidth*1.5, _singleWidth*4);
let _numColors = rndint(5, 15);
let colorsgrid = chroma.scale([colorsshuffle[0], colorsshuffle[1]]).mode('lch').colors(_numColors);
if(fxrand() > 0.95 ) {
    colorsgrid = chroma.scale([colorsshuffle[0], colorsshuffle[0]]).mode('lch').colors(_numColors);
}
let positiverand = rndint(1, 25);
let _noiseAmplitude = rndint(30, 200);
let _noiseCol = rndint(30, 80)/10000;
let _noisePos = rndint(1, 4)/1000;
let _fontSize = 60;
let border;
let unique_contour_values;
let n_rendered_points;
function preload() {
    /*     Ensure the .ttf or .otf font stored in the assets directory
        is loaded before setup() and draw() are called */
    font = loadFont('../assets/fonts/Despairs-X3Wxo.ttf');
}
function setup(){
    colorMode(HSB, 360, 100, 100, 100)
    if(urlParams.has('density')) {
        let pixetemp = int(urlParams.get('density'))
        if(pixetemp > 0 && pixetemp <= 4) {
        pixeld = pixetemp
        console.log("pixelDensity: " + pixeld)
        }
    }  
    //setup text things
    textFont(font);
    textSize(_fontSize);
    pixelDensity(pixeld);
    angleMode(DEGREES);
    noLoop()
    noStroke()
    rectMode(CENTER);
    //seed noise
    seed = int(fxrand()*9999999);
    noiseSeed(seed);
    //setup coloring and render settings
    drawcanvas = createCanvas(_width, _height);
    drawcanvas.parent("fulllscreen");

    //render beautiful background
    background('#f3f3e7');
    bgc =  chroma("#" + colorsshuffle[0]).darken(1).alpha(0.03).hex();
    fill(bgc);
    rect(_width/2, _height/2, _width, _height);
    //define grid of noise
    border = rndint(150, 250);
    border = 200;
    let _numX = Math.round((_width-border) / _singleWidth)
    let _numY = Math.round((_height-border) / _singleHeight)
    /*  we now need to generate the noise on a 2D grid
        this noise will be used to identify the colors and contours */
    noiseSeed(seed+_noiseCounter)
    let noiseColor = [];
    let posx = [];
    let posy = []; //need these to get the bounding box in the same way
    for (let i=0;i<_numX;i++){
        noiseColor[i] = [];
        for (let j=0; j<_numY; j++){
            noiseColor[i][j] = noise(i*_noiseCol, j*_noiseCol);
            posx.push((i*_singleWidth) + (noise(i*_noisePos, j*_noisePos)*_noiseAmplitude) + rndint(-positiverand, positiverand));
            posy.push((j*_singleHeight) + (noise(i*_noisePos, j*_noisePos)*_noiseAmplitude) + rndint(-positiverand, positiverand));
        }
    }

    /*  now that we have random noise generated on a grid, we see how many
        contours that gives us
        We flatten the 2D array to get the min and max values easily */
    let noiseColor_min = noiseColor.reduce(function (p, c){return p.concat(c);}).min();
    let noiseColor_max = noiseColor.reduce(function (p, c){return p.concat(c);}).max();
    // console.log(noiseColor_min,noiseColor_max);
    //now we rescale the noisevalues to discrete steps to identify the number of contours
    //to place
    let noiseContours = [];
    for (let i=0;i<_numX;i++){
        noiseContours[i] = [];
        for (let j=0; j<_numY; j++){
            noiseContours[i][j] = Math.round(map(noise(i*_noiseCol, j*_noiseCol), _noiseFloor, 1, 1, 100));
            //only allow discrete number of colors for the noise
            noiseColor[i][j] = Math.round(map(noiseColor[i][j],noiseColor_min, noiseColor_max, 0, _numColors-1));
        }
    }
    // console.log(_numColors);
    // let unique_color_indices = noiseColor.reduce(function(p,c){return p.concat(c);}).filter(onlyUnique);
    // console.log(unique_color_indices);
    // we now need to get an estimate the volume of each contour to then adjust
    // the number of points populated such that each contour has roughly uniform density
    // start with enforcing uniform density, then maybe can scale based on volume, etc
    unique_contour_values = noiseContours.reduce(function(p,c){return p.concat(c);}).filter(onlyUnique);
    let contour_volumes = noiseContours.reduce(function(p,c){return p.concat(c);}).reduce(
        function (count, currentValue) {
            return (count[currentValue] ? ++count[currentValue] : (count[currentValue] = 1), count);
        },{}); //dictionary not array!

    /*up to now on a 2D grid, we have generated integer level sets of contours
    and the colors associated with that portion of the grid
    now we need to go through each integer level set, generate random points across the canvas,
    if the point is in the given level set, add it with the nearest color value. Easy, right?*/

    let xstart = (_width-posx.max()-posx.min())/2
    let ystart = (_height-posy.max()-posy.min())/2
    // console.log(unique_contour_values);
    let ideal_density = 3.2*(12/unique_contour_values.length); //test case emma: largest contour is ~530 pixels, so we want ~2000points, giving density of 3.7 point/pixel
    n_rendered_points = [];
    for (idn=0;idn<unique_contour_values.length;idn++){
        console.log('Generating contour ' + idn + '/' + unique_contour_values.length);
        let contour_val = unique_contour_values[idn];
        let accumulator = 0;
        let n_rendered = Math.floor(ideal_density*contour_volumes[contour_val]); //total number of points to render for each contour
        n_rendered_points.push(n_rendered);
        while (accumulator < n_rendered){
            let tx = rndfloat(posx.min(), posx.max());
            let ty = rndfloat(posy.min(), posy.max());
            let itx = int(tx/_singleWidth)-1;
            let ity = int(ty/_singleHeight)-1;
            if (itx>=noiseContours.length){itx=noiseContours.length-1;}else if (itx<0){itx=0;}
            if (ity>=noiseContours[0].length){ity=noiseContours[0].length-1;}else if (ity<0){ity=0;}
            // console.log('tx: '+tx+' ty: ' + ty + ' itx: '+itx+'/'+noiseContours.length+' ity: '+ ity+'/'+noiseContours[0].length);
            // console.log(noiseContours[itx][ity]);
            // console.log(,noiseContours[0].length);
            if (
                noiseContours[itx][ity] == contour_val
                && inBounds(tx,ty,posx,posy)
            ){
                let colorindex = noiseColor[itx][ity];
                pointarray.push(new PointObj(tx+xstart,ty+ystart, contour_val,colorindex));
                accumulator++;
            }
        }
    }
    console.log('Generated points');
//     var randomvariables = {
//         BaseColors:String(colors),
//         SingleDistance:String(singlew),
//         Grid:String(aantalw + " x " + aantalh),
//         SpaceStart:String(spacestart),
//         SpaceBetween:String(spacebetween),
//         PosRandom:String(posrand),
//    }
//    console.table(randomvariables);

//    window.$fxhashFeatures = {
//         "Base Colors": String(colors),
//         "Single Distance":String(singlew),
//         "Grid":String(aantalw + " x " + aantalh),
//         "Space Start":String(spacestart),
//         "Space Between":String(spacebetween),
//         "Pos Random":String(posrand),
//    }
//     /* so now we should have n_rendered*unique_contour_values.length points
//     randomly generated and organized in a 2D array (pointarray), with
//     first dimension iterating over the number of contours, and the second
//     over the randomly generated points. Its elements are objects containing
//     their position on the screen, the contour they belong to, and their color*/
}

function draw(){
    blendMode(MULTIPLY);
    let contour = 0;
    let accumulator = 0;
    let y_min = height;
    let id_y_min = pointarray.length-1;
    for (let idp=0;idp<pointarray.length;idp++){
        if (accumulator == n_rendered_points[contour]){
            accumulator = 0;
            contour++;
            console.log('Rendering contour '+contour+'/'+unique_contour_values.length);
        }
        if (pointarray[idp].vec.y < y_min){y_min = pointarray[idp].vec.y;id_y_min=idp;}
        pointarray[idp].createConnections(pointarray);
        pointarray[idp].display();
        accumulator++;
    }
    push();
    let write_col = chroma(colorsgrid[pointarray[id_y_min].colorindex]).alpha(rndint(5, 8)/10).hex();
    fill(write_col);
    textAlign(CENTER, CENTER);
    text(stringSeed.split("").join(" "), width/2,min(95,pointarray[id_y_min].vec.y));
    pop();
    const loadingdiv = document.getElementById('loadingd')
    loadingdiv.remove();
    callpreviewtimer = setTimeout(callpreview, 1000);
}

function callpreview() {
    clearTimeout(callpreviewtimer);
    console.log("done")
    fxpreview()
}
function onlyUnique(value, index, self) {
    return self.indexOf(value) === index;
}

function inBounds(ex, why, posxs, posys){
    if (ex > posxs.min() && ex < posxs.max() && why > posys.min() && why < posys.max()){
        return true;
    } else {
        return false;
    }
}
class PointObj{
    constructor(x, y, contour, colorindex){
        this.vec = createVector(x,y);
        this.contour = contour;
        this.colorindex = colorindex;
        this.connections = [];
    }
    createConnections(pts){
        for (let k=0;k<pts.length;k++){
            if (fxrand() > _noiseFloor){
                if((this != pts[k]) && this.contour == pts[k].contour){
                    if (fxrand() > _singleWidthDiv ){
                        let dist = this.vec.dist(pts[k].vec);
                        if (dist > _minDist && dist < _minDist+_maxLen && this.vec.y > 1.5*border){
                            this.connections.push(
                                [
                                    rndint(_singleWidth/2,_singleWidth)/40, //strokeWeight
                                    chroma(colorsgrid[this.colorindex]).alpha(rndint(5, 8)/10).hex(),//strokecolor
                                    pts[k],//the connection itself
                                ]
                            )
                        }
                    }
                }
            }
        }
    }
    display(){
        for (let k=0;k<this.connections.length;k++){
            // strokeWeight(this.connections[k][0]);
            stroke(this.connections[k][1])
            strokeWeight(0.1);
            // stroke(0);
            let connection = this.connections[k][2].vec;
            line(this.vec.x,this.vec.y,connection.x,connection.y)
        }
    }
}

function randomarray(inputarray) {
    return inputarray[Math.floor(fxrand()*inputarray.length)];
}
function shufflearr(array) {
    for (let i = array.length - 1; i > 0; i--) {
        const j = Math.floor(fxrand() * (i + 1));
        [array[i], array[j]] = [array[j], array[i]];
    }
    return array;
}

Array.prototype.max = function() {
    return Math.max.apply(null, this);
};

Array.prototype.min = function() {
    return Math.min.apply(null, this);
};

function rndint(min, max) {
    return Math.floor(fxrand() * (max - min + 1) + min)
}

function rndfloat(min, max){
    return fxrand()*(max-min+1) + min;
}

function cyrb128(str) {
    let h1 = 1779033703, h2 = 3144134277,
        h3 = 1013904242, h4 = 2773480762;
    for (let i = 0, k; i < str.length; i++) {
        k = str.charCodeAt(i);
        h1 = h2 ^ Math.imul(h1 ^ k, 597399067);
        h2 = h3 ^ Math.imul(h2 ^ k, 2869860233);
        h3 = h4 ^ Math.imul(h3 ^ k, 951274213);
        h4 = h1 ^ Math.imul(h4 ^ k, 2716044179);
    }
    h1 = Math.imul(h3 ^ (h1 >>> 18), 597399067);
    h2 = Math.imul(h4 ^ (h2 >>> 22), 2869860233);
    h3 = Math.imul(h1 ^ (h3 >>> 17), 951274213);
    h4 = Math.imul(h2 ^ (h4 >>> 19), 2716044179);
    return [(h1^h2^h3^h4)>>>0, (h2^h1)>>>0, (h3^h1)>>>0, (h4^h1)>>>0];
}

function exportPNG() {
    saveCanvas(stringSeed, "png");
}
function keyReleased() {
    if (key == 's' || key == 'S' || key == 'e' || key == 'E') exportPNG();
}