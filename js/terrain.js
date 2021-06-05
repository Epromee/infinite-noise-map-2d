
//==========================================================================================

/* Deterministic 2d random - bicycle edition */
function det_random_2d(seed, layers, layer, x, y) {

    // merge coords with layers
    if (x % 2 == y % 2) {
        x = x * layers;
        y = y * layers;
        x = x + layer;
    }
    else {
        x = x * layers;
        y = y * layers;
        y = y + layer;
    }
    
    //-------------------------------------------------------------------------
    
    // compute
    let val = 0;

    // temporary vars for linear transformations
    let nx, ny;

    x += seed;
    
    // xor shuffle
    y ^= x;
    nx = x - y;
    ny = x + y - seed;
    x = nx;
    y = ny;
    x ^= y;

    // radius shuffle
    y ^= ((x * x) & -Math.abs(x)) * Math.sign(x);
    x ^= ((y * y) & -Math.abs(y)) * Math.sign(y);
    y ^= ((x * x) & -Math.abs(x)) * Math.sign(x);
    x ^= ((y * y) & -Math.abs(y)) * Math.sign(y);
    
    val = x ^ y ^ seed;
   
    //-----------------------------------------------------------
    
    // and return
    if (Math.floor(val / 997) % 2 == 0)
        val = ((997 + (val % 997)) % 997) / 997;
    else
        val = ((997 + (-val % 997)) % 997) / 997;
       
    return val;
}

//==========================================================================================

const HALF_PI = Math.PI / 2;
const TWO_PI = Math.PI * 2;

function smoothf_norm(val) {
    return Math.sin(val * HALF_PI);
}

// smooth
function smoothf(val) {

    let fixval = val * 2 - 1;
    let smoothval = fixval;

    smoothval = smoothf_norm(smoothval);

    let backval = (smoothval + 1) / 2;

    return backval;
};

let pnoise_cache = [];
let pnoise_cache_entries = 1000;

/* Let us not recompute the values again each time */
for (let k = 0; k < pnoise_cache_entries; ++k) {
    pnoise_cache.push({ ci : 0, cj : 0, cmseed : 0, crsize : 0, value : 0});
}

/* Perlin-like noise, but using a custom random function (det_random_2d) */
function pnoise(i, j, mseed, rsize) {

    // check cache
    let hash = Math.abs(((rsize << 8) ^ (mseed << 4) ^ (i << 2) ^ (j) ^ (-1))) % pnoise_cache_entries;

    let def = pnoise_cache[hash];

    // if cached, then return the value immediately
    if (def.crsize === rsize && def.cmseed === mseed && def.ci === i && def.cj === j) {
        return def.value;
    }
    
    rsize = Math.floor(rsize);
    
    // if the radius size is too small for any particular randomization
    if (rsize === 0)
        return 0.5;
    
    i = Math.floor(i);
    j = Math.floor(j);

    let idrsize = Math.floor(i / rsize);
    let jdrsize = Math.floor(j / rsize);

    // properties of perlin vectors
    let tl_x, tl_y, bl_x, bl_y, tr_x, tr_y, br_x, br_y;

    // compute angles of the perlin vectors
    let tl_r = det_random_2d(mseed, 1, 0, idrsize, jdrsize);
    let tr_r = det_random_2d(mseed, 1, 0, idrsize + 1, jdrsize);
    let bl_r = det_random_2d(mseed, 1, 0, idrsize, jdrsize + 1);
    let br_r = det_random_2d(mseed, 1, 0, idrsize + 1, jdrsize + 1);
            
    // compute projections of the perlin vectors
    tl_x = Math.sin(TWO_PI * tl_r);
    tr_x = Math.sin(TWO_PI * tr_r);
    bl_x = Math.sin(TWO_PI * bl_r);
    br_x = Math.sin(TWO_PI * br_r);

    tl_y = Math.cos(TWO_PI * tl_r);
    tr_y = Math.cos(TWO_PI * tr_r);
    bl_y = Math.cos(TWO_PI * bl_r);
    br_y = Math.cos(TWO_PI * br_r);
    
    let ceili = Math.ceil(i / rsize);

    // if i is zero, then we can't ceil to the rightmost side, so we do it manually
    if (ceili == idrsize)
        ceili++;

    let ceilj = Math.ceil(j / rsize);
    
    // if i is zero, then we can't ceil to the lowest side, so we do it manually
    if (ceilj == jdrsize)
        ceilj++;

    let idr_r_size = idrsize * rsize;
    let jdr_r_size = jdrsize * rsize;

    let dx = (i - idr_r_size) / (ceili * rsize - idr_r_size);
    let dy = (j - jdr_r_size) / (ceilj * rsize - jdr_r_size);

    let dxm = dx - 1;
    let dym = dy - 1;

    // select all four powers to the perlin vectors depending on the dot product
    let tl = (dx * tl_x) + (dy * tl_y);
    let tr = (dxm * tr_x) + (dy * tr_y);
    let bl = (dx * bl_x) + (dym * bl_y);
    let br = (dxm * br_x) + (dym * br_y);

    dx = smoothf(dx);
    dy = smoothf(dy);

    tl = smoothf_norm(tl);
    tr = smoothf_norm(tr);
    bl = smoothf_norm(bl);
    br = smoothf_norm(br);

    dxm = 1 - dx;
    dym = 1 - dy;

    // and then select one particular approximated power, depending on the distances
    let power = (tl * dxm + tr * dx) * dym + (bl * dxm + br * dx) * dy;

    // map [-1; 1] to [0; 1]
    let val = (power + 1) / 2;

    // write cache back
    def.ci = i;
    def.cj = j;
    def.cmseed = mseed;
    def.crsize = rsize;
    def.value = val;
    
    return val;
}

/* This function draws a certain terrain pixel */
function drawBeautifulPixel(i, j) {

    /* [DISCLAIMER!] Sorry for the magical numbers, I had to select the values from nowhere using only my experiments */

    let seed = 6661337;

    let level_temperature_ratio = 0.05;
    
    // P variables are noise sources which I combine in some ways (no strict way) to achieve the needed effect
    let p0 = pnoise(i, j, seed, 640);
    let p1 = pnoise(i, j, seed, 320);
    let p2 = pnoise(i, j, seed, 160);
    let p3 = pnoise(i, j, seed, 80);
    let p4 = pnoise(i, j, seed, 40);
    let p5 = pnoise(i, j, seed, 20);
    let p6 = pnoise(i, j, seed, 10);
    let p7 = pnoise(i, j, seed, 5);
    let p8 = pnoise(i, j, seed, 2);
    
    let q1 = 0, q2 = 0, q3 = 0, q4 = 0, q5 = 0, q6 = 0;
    
    // well this is kinda "height" level
    let lvl = (p0 / 1 + p1 / 2 + p2 / 3 + p3 / 4 + p4 / 5 + p5 / 6 + p6 / 7 + p7 / 8 + p8 / 9) / (1 + 1 / 2 + 1 / 3 + 1 / 4 + 1 / 5 + 1 / 6 + 1 / 7 + 1 / 8 + 1 / 9);
    
    // this is "temperature" level
    let temp = 0;

    let r = 0, g = 0, b = 0;

    // Q variables are just like P variables in case I need more noise sources
    let calc_q = function () {
        q1 = pnoise(i, j, seed + 2, 300);
        q2 = pnoise(i, j, seed + 3, 88);
        q3 = pnoise(i, j, seed + 5, 44);
        q4 = pnoise(i, j, seed + 7, 14);
        q5 = pnoise(i, j, seed + 11, 6);
        q6 = pnoise(i, j, seed + 13, 2);
        
        // well, these guys recompute the temperature depending on the temperature itself and the destinations to the polar sides of the map
        let map_height = (j / canH);
        let noise_ratio = (q1 / 2 + q2 / 3 + q3 / 4 + q4 / 5 + q5 / 5 + q6 / 5) / (1 / 2 + 1 / 3 + 1 / 4 + 1 / 5 + 1 / 5 + 1 / 5);

        let pr = 4 * map_height * map_height - 4 * map_height + 1; // [1 ... 0 ... 1] function
        let prn = 1 - pr;   // [0 ... 1 ... 0] function

        let rs = Math.pow((pr * prn), 0.15);

        temp = (1 - rs) * prn + rs * noise_ratio;
    };

    // ocean
    if (lvl < 0.5) {
        
        // water
        b = 127 * lvl;

        // underwater volcanos
        let can0 = pnoise(i, j, seed + 666, 20) * 2 - 1;
        let can1 = pnoise(i, j, seed + 666, 15) * 2 - 1;

        if (Math.abs(can0 * can1) < (0.5 - lvl * 1.2)) {
            
            b = b * 0.8 + 20 * lvl * 0.2
            r = g = 5 * lvl;

            if (Math.abs(can0 * can1) < (0.5 - lvl * 1.52)) {
                r = 20;
            }
        }
    }
    
    // grass / sand
    else if (lvl >= 0.5 && lvl < 0.6) {

        calc_q();

        let inr_lvl = (lvl - 0.5) / (0.7 - 0.5);
        //temp = Math.min(temp, inr_lvl * 0.45 + 0.55); 
        temp = Math.min(temp, 1 - inr_lvl * inr_lvl);
        temp = Math.min(temp, 0.55 + inr_lvl);        

        let river_drawn = false;

        // rivers
        if (temp >= 0.4) {
            let riv0 = pnoise(i, j, seed + 777, 80) * 2 - 1;
            let riv1 = pnoise(i, j, seed + 777, 20) * 2 - 1;
            let riv2 = ((lvl - 0.5) / (lvl - 0.6)) * 2 - 1;

            let riv_lvl = Math.abs(riv0 * riv1 * riv2);
            let riv_threshold = (0.56 - lvl) * 0.15;

            if (riv_lvl < riv_threshold) {

                let approx_gate = 0.55;
                let lvlr = (lvl - 0.5) / (approx_gate - 0.5);

                r = (lvl < approx_gate ? (10 * lvlr + 0 * (1 - lvlr)) : 10);
                g = (lvl < approx_gate ? (60 * lvlr + 0 * (1 - lvlr))  : 60);
                b = (lvl < approx_gate ? (160 * lvlr + 127 / 2 * (1 - lvlr)) : 160);
                river_drawn = true;
            }

            temp = Math.min(temp, (riv_lvl - riv_threshold) * 10 + 0.57);
        }

        // if no river, then
        if (!river_drawn) {
            
            // snow
            if (temp < 0.4) {
                r = g = b = 255;
                if (temp > 0.36) {
                    g = 20 + 60 * (1 - (temp - 0.36) / (0.4 - 0.36));
                    r = 30 * (1 - (temp - 0.36) / (0.4 - 0.36));
                    b = r;
                }
            }
            // grass
            else if (temp >= 0.4 && temp < 0.6) {
                g = 20 + 80 * (temp - 0.4) / (0.6 - 0.4);
            }
            //sand
            else if (temp > 0.6) {
                r = 255 * 4 * (temp - 0.6) / (1 - 0.6);
                g = 80 + 100 * 4 * (temp - 0.6) / (1 - 0.6);
                r = Math.max(0, Math.min(r, 200));
                g = Math.max(0, Math.min(g, 200));
                b = Math.min(g / 2, r / 2);
            }
        }
    }
       
    // mountain
    else if (lvl >= 0.6) {
        
        calc_q();
        
        let inr_lvl = (lvl - 0.5) / (0.7 - 0.5);
        temp = Math.min(temp, 1 - inr_lvl * inr_lvl);   

        let rock_look = (p7 + p8 + q5 * 2 + q6 * 3) / 7;
        
        r = g = b = 50 + (180 * Math.sqrt((lvl - 0.6)) / (1 - 0.6) ) * rock_look;

        // mountain snow
        if (temp < 0.4 && lvl < 0.71) {
            r = g = b = 255;
        }

        // mountain volcano
        if (lvl > 0.73) {
            r = 255;
            g = 120;
            b = 0;
        }
    }

    r = Math.max(0, Math.min(r, 255));
    g = Math.max(0, Math.min(g, 255));
    b = Math.max(0, Math.min(b, 255));

    return { r : r, g : g, b : b};
}

var canvas = document.getElementById("myCanvas");

var ctx = canvas.getContext("2d");

let canW = canvas.clientWidth;
let canH = canvas.clientHeight;

ctx.rect(0, 0, canW, canH);
ctx.fill();

let imageData = ctx.getImageData(canW - 1, 0, 1, canH);
let imageDataData = imageData.data;

let iterations = 0;

function redraw() {

    iterations++;
    
    for (let j = 0; j < canH; j += 1) {
        let rgb = drawBeautifulPixel(iterations, j);
        imageDataData[((j) << 2) + 0] = rgb.r;
        imageDataData[((j) << 2) + 1] = rgb.g;
        imageDataData[((j) << 2) + 2] = rgb.b;
    };
    
}

function sec() {
    ctx.putImageData(ctx.getImageData(1, 0, canW - 1, canH), 0, 0);
    redraw();
    ctx.putImageData(imageData, canW - 1, 0);
    window.requestAnimationFrame(sec);
}

window.requestAnimationFrame(sec);
