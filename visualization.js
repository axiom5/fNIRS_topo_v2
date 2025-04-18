// visualization.js

// Ensure EEG_DATA is defined (should be injected by Python before this script runs)
if (typeof EEG_DATA === 'undefined') {
    console.error("CRITICAL: EEG_DATA is not defined. Data injection likely failed.");
    // Optionally display an error message to the user
    document.body.innerHTML = '<h1>Error: Failed to load visualization data.</h1>';
}

const fnirsVisApp = {
    // --- Configuration & Constants ---
    config: {
        GRID_SIZE: 80,
        MIN_WEIGHT_SUM_THRESHOLD: 1e-3,
        UNKNOWN_REGION_NAME: "Unknown", // Consistent with Python
        OUTSIDE_BRAIN_NAME: "Outside",   // Consistent with Python
        TIMESERIES_HEIGHT: 150,
        TIMESERIES_MARGIN: { top: 10, right: 20, bottom: 30, left: 50 },
        HEAD_CLIP_ID: "head-clip",
        DEFAULT_SUBJECT_RANGE: [-1, 1],
        // Add CSS selectors if needed, e.g., TOPO_SVG_SELECTOR: '#topomap'
    },

    REGION_COLORS: {
        // Finer PFC (Blue/Teal Spectrum)
        "Frontopolar": "#5966a6",             // Deeper Blue (Frontopolar)
        "DLPFC": "#1f78b4",          // Medium Blue (Dorsolateral)
        // -- Revised VLPFC/IFG --
        "VLPFC/IFG": "#4db4ab",      // Teal/Cyanish-Blue

        "Dorsal PFC": "#b2df8a",      // Light Green (Keeping distinct, spatially dorsal)
        "OFC": "#fdbe85",            // Light Orange/Peach (Orbitofrontal)
        "VMPFC": "#bae4bc",          // Light Mint Green (Ventromedial)

        // Cingulate (Pinks/Purples)
        "ACC": "#e377c2",            // Pink (Anterior)
        "Post Cingulate": "#f7b6d2", // Lighter Pink (Posterior)

        // Motor/Premotor (Reds/Oranges)
        "Premotor": "#e31a1c",          // Bright Red
        "Motor": "#fdbf6f",       // Light Orange/Yellow (Link to Motor, but distinct)

        // Other regions (Distinct Colors)
        "Somatosensory": "#ff7f00",  // Orange (Distinct from Premotor's new color)
        "Parietal": "#b15928",       // Brown
        "Temporal": "#6a3d9a",       // Dark Purple
        "Temporal (Aud)": "#cab2d6", // Light Lilac
        "Temporal (Ins)": "#9e9ac8", // Grayish Purple (Insular) - Adjusted for distinction
        "Occipital": "#ffff99",      // Light Yellow
        "Limbic": "#7f7f7f",          // Gray

        // Fallbacks
        "Unknown": "#bdbdbd",        // Medium Gray
        "Outside": "#f0f0f0",        // Very Light Gray
    },

    // --- State Variables ---
    state: {
        currentSubjectId: null,
        currentSubjectData: null, // { times: [], data: [[]] }
        currentSubjectMin: -1,
        currentSubjectMax: 1,
        nTimes: 0,
        nAllChannels: 0,
        currentActiveChannelIndices: [], // Indices relative to EEG_DATA.channels.names
        channelMidpointSvgPos: [], // Array of {x, y, channelIndex, channelName}
        allChannelRegions: [], // Stores the ba_regions array from Python
        currentVisualizationMode: 'heatmap', // 'heatmap' or 'anatomical'
        selectedChannelIndex: null,
        selectedChannelName: null,
        gridPoints3D: [], // Precomputed [x,y,z] in Digitizer space for interpolation
        // Geometry calculated during init
        svgWidth: 0,
        svgHeight: 0,
        centerX: 0,
        centerY: 0,
        headRadiusPixels: 0,
        xMinGrid: 0, xMaxGrid: 0, yMinGrid: 0, yMaxGrid: 0,
        gridStepX: 0, gridStepY: 0,
        tsWidth: 0, tsHeight: 0,
    },

    // --- D3 Selections & Scales ---
    elements: {
        // Controls
        timeSlider: null, timeDisplay: null,
        smoothingSlider: null, smoothingDisplay: null,
        subjectSelect: null, subgroupSelect: null, viewModeSelect: null,
        timeseriesLabel: null,
        // SVG Main Groups
        svg: null, gMain: null,
        heatmapGroup: null, channelLinesGroup: null, optodeDotsGroup: null,
        channelClickTargetsGroup: null, fiducialsGroup: null, anatomicalRegionsGroup: null,
        // Timeseries Plot
        timeseriesSvg: null, tsXAxisGroup: null, tsYAxisGroup: null, tsLinePath: null, tsCurrentTimeMarker: null,
        // Colorbar
        colorbarSvg: null,
        // Legend
        legendContainer: null,
    },
    scales: {
        colorScaleSignal: null, // D3 scale for signal heatmap
        regionColorScale: null, // D3 scale for regions
        tsXScale: null, tsYScale: null,
    },

    // --- Initialization ---
    init: function() {
        console.log("Initializing fNIRS Visualization App...");
        if (!this._validateInjectedData()) {
            this._showFatalError("fNIRS data (EEG_DATA) is invalid or missing essential fields. Check console.");
            return;
        }
        this._updateTitles();
        this._cacheDomElements();
        this._setupSvg();
        this._calculateGeometry(); // Depends on SVG setup
        this._setupScales();
        this._drawStaticHeadElements();
        this._calculateSvgPositions(); // Depends on geometry
        this._drawOptodesAndChannels(); // Depends on SVG positions
        this._drawClickTargets(); // Depends on SVG positions
        this._drawFiducials(); // Depends on geometry
        this._precomputeGrid3D(); // Depends on geometry
        this._setupColorbar();
        this._setupTimeseriesPlot();
        this._setupViewModeSelection(); // Setup controls before loading data
        this._setupSubgroupSelection();
        this._setupSubjectSelection(); // Loads first subject, triggers initial updates
        this._bindEventListeners();

        console.log("Initialization complete.");
    },

    _validateInjectedData: function() {
        const d = EEG_DATA; // Alias for brevity
        if (!d || Object.keys(d).length === 0) { console.error("EEG_DATA is empty."); return false; }
        if (!d.subjects || !d.subjectList) { console.error("Missing subjects/subjectList."); return false; }
        if (!d.geometry) { console.error("Missing geometry."); return false; }
        if (!d.channels || !d.channels.names || !d.channels.positions3d || !d.channels.ba_regions) { console.error("Missing channel info (names, positions3d, ba_regions)."); return false; }
        if (!Array.isArray(d.channels.ba_regions)) { console.error("channel.ba_regions is not an array."); return false; }
        if (!d.optodes || !d.optodes.pairs) { console.error("Missing optodes info (pairs)."); return false; }
        // Basic length check (more thorough checks happen during processing)
        if (d.channels.names.length !== d.channels.positions3d.length || d.channels.names.length !== d.channels.ba_regions.length) {
             console.error("Channel data arrays have inconsistent lengths."); return false;
        }
        if (d.channelSubgroups === undefined || d.defaultSubgroup === undefined) { console.error("Missing channelSubgroups or defaultSubgroup."); return false; }
        // Add more checks as needed (e.g., check specific geometry fields)
        return true;
    },

    _showFatalError: function(message) {
        console.error("FATAL ERROR:", message);
        const container = document.querySelector('.container');
        if (container) {
            container.innerHTML = `<h1>Error</h1><p>${message}</p>`;
        } else {
            document.body.innerHTML = `<h1>Error</h1><p>${message}</p>`;
        }
    },

    _updateTitles: function() {
        document.title = `fNIRS Topo: T${EEG_DATA.trialNum} ${EEG_DATA.chromophore}`;
        document.getElementById('main-title').textContent = `fNIRS Topography: Trial ${EEG_DATA.trialNum}`;
        document.getElementById('subtitle').textContent = `${EEG_DATA.chromophore} (${EEG_DATA.channelType})`;
    },

    _cacheDomElements: function() {
        this.elements.timeSlider = document.getElementById('time-slider');
        this.elements.timeDisplay = document.getElementById('time-display');
        this.elements.smoothingSlider = document.getElementById('smoothing-slider');
        this.elements.smoothingDisplay = document.getElementById('smoothing-display');
        this.elements.timeseriesLabel = document.getElementById('timeseries-label');
        this.elements.subjectSelect = document.getElementById('subject-select');
        this.elements.subgroupSelect = document.getElementById('subgroup-select');
        this.elements.viewModeSelect = document.getElementById('view-mode-select');
        this.elements.legendContainer = document.getElementById('anatomical-legend');
        this.elements.colorbarSvg = d3.select('#colorbar');
    },

    _setupSvg: function() {
        const topoContainer = d3.select('.topomap-container');
        if (topoContainer.empty()) { console.error("Cannot find .topomap-container"); return; }
        const containerWidth = topoContainer.node().getBoundingClientRect().width;
        this.state.svgWidth = containerWidth;
        this.state.svgHeight = this.state.svgWidth * (500/600); // Maintain aspect ratio

        this.elements.svg = d3.select('#topomap')
            .attr('width', this.state.svgWidth)
            .attr('height', this.state.svgHeight)
            .attr('viewBox', `0 0 ${this.state.svgWidth} ${this.state.svgHeight}`)
            .attr('preserveAspectRatio', 'xMidYMid meet');
        this.elements.svg.selectAll('*').remove(); // Clear previous content

        const margin = { top: 50, right: 50, bottom: 50, left: 50 };
        this.elements.gMain = this.elements.svg.append('g').attr('transform', `translate(${margin.left}, ${margin.top})`);

        // Define clip path
        this.elements.gMain.append("defs").append("clipPath").attr("id", this.config.HEAD_CLIP_ID)
           .append("circle")
           .attr("id", "clip-circle"); // Give circle an ID for geometry calculation

        // Add drawing layers
        this.elements.heatmapGroup = this.elements.gMain.append("g").attr("class", "heatmap-group").attr("clip-path", `url(#${this.config.HEAD_CLIP_ID})`);
        this.elements.anatomicalRegionsGroup = this.elements.gMain.append("g").attr("class", "anatomical-regions-group").attr("clip-path", `url(#${this.config.HEAD_CLIP_ID})`).style("display", "none"); // Hidden initially
        this.elements.channelLinesGroup = this.elements.gMain.append("g").attr("class", "channel-lines-group");
        this.elements.optodeDotsGroup = this.elements.gMain.append("g").attr("class", "optode-dots-group");
        this.elements.fiducialsGroup = this.elements.gMain.append("g").attr("class", "fiducials-group");
        this.elements.channelClickTargetsGroup = this.elements.gMain.append("g").attr("class", "channel-click-targets");
    },

    _calculateGeometry: function() {
        const clipCircle = this.elements.gMain.select("#clip-circle");
        if (clipCircle.empty()) { console.error("Clip circle not found for geometry calculation."); return; }

        const margin = { top: 50, right: 50, bottom: 50, left: 50 }; // Re-access margin needed? Define higher up?
        const innerWidth = this.state.svgWidth - margin.left - margin.right;
        const innerHeight = this.state.svgHeight - margin.top - margin.bottom;
        this.state.headRadiusPixels = Math.min(innerWidth, innerHeight) / 2 * 0.95;
        this.state.centerX = innerWidth / 2;
        this.state.centerY = innerHeight / 2;

        // Update clip path circle dimensions now that geometry is known
        clipCircle.attr("cx", this.state.centerX).attr("cy", this.state.centerY).attr("r", this.state.headRadiusPixels);
    },

    _setupScales: function() {
        this.scales.colorScaleSignal = d3.scaleSequential(d3.interpolateTurbo);
        this.scales.regionColorScale = d3.scaleOrdinal(); // Domain set in setupLegend
        this.scales.tsXScale = d3.scaleLinear();
        this.scales.tsYScale = d3.scaleLinear();
    },

    _drawStaticHeadElements: function() {
        const { centerX, centerY, headRadiusPixels } = this.state;
        // Head Outline
        this.elements.gMain.append('circle')
            .attr('class', 'head-outline')
            .attr('cx', centerX).attr('cy', centerY).attr('r', headRadiusPixels);
        // Nose
        if (EEG_DATA.head && EEG_DATA.head.nose) {
             const { x: nose_norm_x, y: nose_norm_y } = EEG_DATA.head.nose;
             const noseAngle = Math.atan2(nose_norm_y, nose_norm_x);
             const noseSize = headRadiusPixels * 0.1;
             const noseWidth = noseSize * 0.5;
             // Calculate base and tip in SVG coords (Y is inverted in SVG)
             const noseBaseSvgX = centerX + headRadiusPixels * Math.cos(noseAngle);
             const noseBaseSvgY = centerY - headRadiusPixels * Math.sin(noseAngle); // SVG Y-axis inversion
             const noseTipSvgX = noseBaseSvgX + noseSize * Math.cos(noseAngle);
             const noseTipSvgY = noseBaseSvgY - noseSize * Math.sin(noseAngle); // SVG Y-axis inversion
             const perpAngle = noseAngle + Math.PI / 2;
             // Calculate points for the triangle path
             const p1x = noseBaseSvgX + noseWidth * Math.cos(perpAngle);
             const p1y = noseBaseSvgY - noseWidth * Math.sin(perpAngle); // SVG Y-axis inversion
             const p2x = noseTipSvgX;
             const p2y = noseTipSvgY;
             const p3x = noseBaseSvgX - noseWidth * Math.cos(perpAngle);
             const p3y = noseBaseSvgY + noseWidth * Math.sin(perpAngle); // SVG Y-axis inversion
             this.elements.gMain.append('path')
                 .attr('class', 'nose-outline')
                 .attr('d', `M ${p1x} ${p1y} L ${p2x} ${p2y} L ${p3x} ${p3y} Z`);
        } else {
            console.warn("Nose data not found in EEG_DATA.head.nose");
        }
    },

    // Helper to calculate SVG coords from normalized projection coords
    _normalizedToSvg: function(norm_x, norm_y) {
        const { centerX, centerY, headRadiusPixels } = this.state;
        return {
            x: centerX + norm_x * headRadiusPixels,
            y: centerY - norm_y * headRadiusPixels // Invert Y for SVG
        };
    },

     _calculateSvgPositions: function() {
        const { headCenter3d, projectionScaleFactor } = EEG_DATA.geometry;
        const channelPoints3D_digi = EEG_DATA.channels.positions3d;
        this.state.channelMidpointSvgPos = []; // Reset

        if (!channelPoints3D_digi || !headCenter3d || !projectionScaleFactor) {
             console.error("Missing data required for channel midpoint SVG projection."); return;
        }

        channelPoints3D_digi.forEach((p3d, index) => {
            if (!p3d || p3d.some(isNaN)) {
                 console.warn(`Skipping channel index ${index} due to invalid position3d.`);
                 // Add placeholder? Maybe just skip. Depends on downstream needs.
                 // For now, skip adding it to state.channelMidpointSvgPos
                 return;
            }
            const proj_x_y = this._projectPointAzimuthalEquidistantJS(p3d, headCenter3d);
            const norm_x = proj_x_y[0] / projectionScaleFactor;
            const norm_y = proj_x_y[1] / projectionScaleFactor;
            const svgPos = this._normalizedToSvg(norm_x, norm_y);
            const channelName = EEG_DATA.channels.names[index];
            this.state.channelMidpointSvgPos.push({ x: svgPos.x, y: svgPos.y, channelIndex: index, channelName: channelName });
        });
    },

    _drawOptodesAndChannels: function() {
        const optodeSvgPos = {}; // Temp map { 'S1': {x, y}, ... }
        // Calculate SVG positions for sources
        Object.entries(EEG_DATA.optodes.sources || {}).forEach(([n, p]) => {
            const svgPos = this._normalizedToSvg(p.x, p.y); optodeSvgPos[n] = svgPos;
        });
        // Calculate SVG positions for detectors
        Object.entries(EEG_DATA.optodes.detectors || {}).forEach(([n, p]) => {
             const svgPos = this._normalizedToSvg(p.x, p.y); optodeSvgPos[n] = svgPos;
        });

        // Draw Channel Lines
        const channelNames = EEG_DATA.channels.names;
        const channelPairs = EEG_DATA.optodes.pairs;
        this.elements.channelLinesGroup.selectAll('.channel-line').remove();
        if (channelNames && channelPairs && channelNames.length === channelPairs.length) {
             channelNames.forEach((name, index) => {
                 const pair = channelPairs[index]; const sName = pair[0], dName = pair[1];
                 if (optodeSvgPos[sName] && optodeSvgPos[dName]) {
                     this.elements.channelLinesGroup.append("line")
                         .datum(index) // Bind original channel index
                         .attr("class", "channel-line")
                         .attr("x1", optodeSvgPos[sName].x).attr("y1", optodeSvgPos[sName].y)
                         .attr("x2", optodeSvgPos[dName].x).attr("y2", optodeSvgPos[dName].y);
                 } else { console.warn(`Missing optode SVG position for channel ${name} (S:${sName}, D:${dName})`); }
             });
        } else { console.error("Channel names/pairs mismatch or missing."); }

        // Draw Optode Dots
        this.elements.optodeDotsGroup.selectAll('.optode').remove();
        Object.entries(optodeSvgPos).forEach(([name, pos]) => {
             const isSrc = name.startsWith('S');
             this.elements.optodeDotsGroup.append("circle")
                 .attr("class",`optode ${isSrc?'source':'detector'}`)
                 .attr("cx",pos.x).attr("cy",pos.y).attr("r", 3);
        });
    },

     _drawClickTargets: function() {
        this.elements.channelClickTargetsGroup.selectAll('.channel-click-target').remove();
        this.state.channelMidpointSvgPos.forEach(midpointData => {
            this.elements.channelClickTargetsGroup.append("circle")
                .datum(midpointData) // Bind { x, y, channelIndex, channelName }
                .attr("class", "channel-click-target")
                .attr("cx", midpointData.x)
                .attr("cy", midpointData.y)
                .on("click", (event, d) => this.handleChannelClick(event, d)); // Use arrow fn
        });
    },

    _drawFiducials: function() {
        this.elements.fiducialsGroup.selectAll('*').remove();
        if (EEG_DATA.fiducials && Object.keys(EEG_DATA.fiducials).length > 0) {
            Object.entries(EEG_DATA.fiducials).forEach(([name, normPos]) => {
                const svgPos = this._normalizedToSvg(normPos.x, normPos.y);
                const markerSize = 4;
                // Draw cross marker
                this.elements.fiducialsGroup.append("path")
                    .attr("class", "fiducial-marker")
                    .attr("d", `M ${svgPos.x - markerSize} ${svgPos.y} L ${svgPos.x + markerSize} ${svgPos.y} M ${svgPos.x} ${svgPos.y - markerSize} L ${svgPos.x} ${svgPos.y + markerSize}`)
                    .style("stroke", "purple").style("stroke-width", 1.5).style("pointer-events", "none");
                // Draw label
                this.elements.fiducialsGroup.append("text")
                    .attr("class", "fiducial-label")
                    .attr("x", svgPos.x + markerSize + 2)
                    .attr("y", svgPos.y + markerSize / 2) // Adjust vertical alignment slightly
                    .style("font-size", "9px").style("fill", "purple").style("pointer-events", "none")
                    .text(name);
             });
        }
    },

    _precomputeGrid3D: function() {
        console.time("Precompute 3D Grid");
        const { centerX, centerY, headRadiusPixels } = this.state;
        const { headRadius3d, headCenter3d, projectionScaleFactor } = EEG_DATA.geometry;
        const GRID_SIZE = this.config.GRID_SIZE;

        this.state.xMinGrid = centerX - headRadiusPixels; this.state.xMaxGrid = centerX + headRadiusPixels;
        this.state.yMinGrid = centerY - headRadiusPixels; this.state.yMaxGrid = centerY + headRadiusPixels;
        this.state.gridStepX = (this.state.xMaxGrid - this.state.xMinGrid) / GRID_SIZE;
        this.state.gridStepY = (this.state.yMaxGrid - this.state.yMinGrid) / GRID_SIZE;
        this.state.gridPoints3D = [];

        const radius3d = headRadius3d > 0 ? headRadius3d : 1.0; // Use DIGITIZER radius

        for (let i = 0; i <= GRID_SIZE; i++) {
            const gridRow3D = [];
            const svg_x = this.state.xMinGrid + i * this.state.gridStepX;
            for (let j = 0; j <= GRID_SIZE; j++) {
                const svg_y = this.state.yMinGrid + j * this.state.gridStepY;
                // Convert SVG back to normalized projection coords
                const norm_x = (svg_x - centerX) / headRadiusPixels;
                const norm_y = -(svg_y - centerY) / headRadiusPixels; // Y inversion
                // Convert normalized proj coords back to angles (Azimuthal Equidistant)
                const proj_x = norm_x * projectionScaleFactor;
                const proj_y = norm_y * projectionScaleFactor;
                const theta = Math.sqrt(proj_x*proj_x + proj_y*proj_y); // Angle from Z axis
                const phi = Math.atan2(proj_y, proj_x); // Angle in XY plane

                let point3d = [NaN, NaN, NaN];
                // Check if point is within the valid angular range for projection
                if (theta <= Math.PI + 1e-6) {
                     // Inverse projection: Convert angles (theta, phi) back to 3D Cartesian offset
                     const dx = radius3d * Math.sin(theta) * Math.cos(phi);
                     const dy = radius3d * Math.sin(theta) * Math.sin(phi);
                     const dz = radius3d * Math.cos(theta);
                     // Add offset to the head center (Digitizer space)
                     point3d = [ headCenter3d[0] + dx, headCenter3d[1] + dy, headCenter3d[2] + dz ];
                }
                gridRow3D.push(point3d);
            }
            this.state.gridPoints3D.push(gridRow3D);
        }
        console.timeEnd("Precompute 3D Grid");
        console.log(`Precalculated 3D grid points (${GRID_SIZE+1}x${GRID_SIZE+1}) for interpolation.`);
    },

    _setupColorbar: function() {
        this.elements.colorbarSvg.selectAll('*').remove(); // Clear previous
        const cbWidth = parseFloat(this.elements.colorbarSvg.style('width')) || 500;
        const cbHeight = parseFloat(this.elements.colorbarSvg.style('height')) || 30;
        const cbInnerWidth = cbWidth * 0.9;
        const cbMargin = { left: cbWidth * 0.05, right: cbWidth * 0.05 };
        const cbDefs = this.elements.colorbarSvg.append('defs');
        const cbGrad = cbDefs.append('linearGradient').attr('id', 'colorbar-gradient').attr('x1', '0%').attr('y1', '0%').attr('x2', '100%').attr('y2', '0%');
        const numStops = 100;
        for (let i = 0; i < numStops; i++) {
            const offset = i / (numStops - 1);
            cbGrad.append('stop').attr('offset', `${offset*100}%`).attr('stop-color', d3.interpolateTurbo(offset));
        }
        this.elements.colorbarSvg.append('rect')
            .attr('x', cbMargin.left).attr('y', 0)
            .attr('width', cbInnerWidth).attr('height', cbHeight * 0.6) // Adjust height of bar
            .style('fill', 'url(#colorbar-gradient)');
    },

    _setupTimeseriesPlot: function() {
        const tsContainerWidth = d3.select('.timeseries-container').node().getBoundingClientRect().width;
        const margin = this.config.TIMESERIES_MARGIN;
        this.state.tsWidth = Math.max(100, tsContainerWidth - margin.left - margin.right);
        this.state.tsHeight = this.config.TIMESERIES_HEIGHT - margin.top - margin.bottom;

        const tsPlotSvg = d3.select("#timeseries-plot")
             .attr('viewBox', `0 0 ${tsContainerWidth} ${this.config.TIMESERIES_HEIGHT}`);

        // Select or create the main group
        this.elements.timeseriesSvg = tsPlotSvg.select("g.main-ts-group");
        if (this.elements.timeseriesSvg.empty()) {
             this.elements.timeseriesSvg = tsPlotSvg.append("g")
                 .attr("class", "main-ts-group") // Add class for selection
                 .attr("transform", `translate(${margin.left},${margin.top})`);

             this.elements.tsXAxisGroup = this.elements.timeseriesSvg.append("g")
                .attr("class", "axis x-axis")
                .attr("transform", `translate(0,${this.state.tsHeight})`);
             this.elements.tsYAxisGroup = this.elements.timeseriesSvg.append("g")
                .attr("class", "axis y-axis");
             // Add X axis label
             this.elements.timeseriesSvg.append("text")
                .attr("class", "x-axis-label")
                .attr("text-anchor", "middle")
                .attr("x", this.state.tsWidth / 2)
                .attr("y", this.state.tsHeight + margin.bottom - 5) // Position below axis
                .style("font-size", "10px").text("Time (s)");
             this.elements.tsLinePath = this.elements.timeseriesSvg.append("path")
                .attr("class", "signal-line");
             this.elements.tsCurrentTimeMarker = this.elements.timeseriesSvg.append("line")
                .attr("class", "current-time-marker").style('visibility', 'hidden');
        } else { // Update existing elements if size changed
             this.elements.tsXAxisGroup.attr("transform", `translate(0,${this.state.tsHeight})`);
             this.elements.timeseriesSvg.select(".x-axis-label").attr("x", this.state.tsWidth / 2);
        }
        // Update scales range
        this.scales.tsXScale.range([0, this.state.tsWidth]);
        this.scales.tsYScale.range([this.state.tsHeight, 0]);
    },

    _setupViewModeSelection: function() {
        if (!this.elements.viewModeSelect) return;
        // Set initial state from HTML default
        this.state.currentVisualizationMode = this.elements.viewModeSelect.value;
        this.elements.viewModeSelect.addEventListener('change', (event) => {
             this.state.currentVisualizationMode = event.target.value;
             console.log("Switched view mode to:", this.state.currentVisualizationMode);
             this._setupLegend(); // Update legend visibility/content
             this.updateTopography(); // Redraw based on new mode
        });
    },

    _setupSubgroupSelection: function() {
        if (!this.elements.subgroupSelect) return;
        this.elements.subgroupSelect.innerHTML = ''; // Clear existing
        // Add "All Channels"
        const allOption = document.createElement('option');
        allOption.value = "All Channels"; allOption.textContent = "All Channels";
        this.elements.subgroupSelect.appendChild(allOption);
        // Add options from EEG_DATA
        if (EEG_DATA.channelSubgroups && Object.keys(EEG_DATA.channelSubgroups).length > 0) {
             Object.keys(EEG_DATA.channelSubgroups).sort().forEach(groupName => {
                 const option = document.createElement('option');
                 option.value = groupName;
                 option.textContent = this._formatGroupName(groupName); // Use helper
                 this.elements.subgroupSelect.appendChild(option);
             });
        } else { console.log("No specific channel subgroups defined in EEG_DATA."); }
        // Set default value
        let defaultVal = "All Channels";
        if (EEG_DATA.defaultSubgroup && (EEG_DATA.defaultSubgroup === "All Channels" || (EEG_DATA.channelSubgroups && EEG_DATA.channelSubgroups[EEG_DATA.defaultSubgroup]))) {
             defaultVal = EEG_DATA.defaultSubgroup;
        } else if (EEG_DATA.defaultSubgroup) { console.warn(`Default subgroup "${EEG_DATA.defaultSubgroup}" not found.`); }
        this.elements.subgroupSelect.value = defaultVal;

        // Listener added in _bindEventListeners
    },

    _setupSubjectSelection: function() {
        if (!this.elements.subjectSelect) return;
        this.elements.subjectSelect.innerHTML = ''; // Clear existing
        if (!EEG_DATA.subjectList || EEG_DATA.subjectList.length === 0) {
             console.error("No subjects found in EEG_DATA.subjectList.");
             // Optionally disable the select or show a message
             return;
        }
        EEG_DATA.subjectList.forEach(subId => {
             const option = document.createElement('option');
             option.value = subId; option.textContent = subId;
             this.elements.subjectSelect.appendChild(option);
        });
        // Load the first subject initially
        this.switchSubject(EEG_DATA.subjectList[0]);

         // Listener added in _bindEventListeners
    },

    _bindEventListeners: function() {
        this.elements.timeSlider?.addEventListener('input', () => this.updateTopography());
        this.elements.smoothingSlider?.addEventListener('input', () => this.updateTopography());
        this.elements.subjectSelect?.addEventListener('change', (event) => this.switchSubject(event.target.value));
        this.elements.subgroupSelect?.addEventListener('change', () => {
            this._updateActiveChannelIndices();
            this.updateTopography(); // Redraw after subgroup change
        });
         // Click listeners for channels are bound in _drawClickTargets
         // View mode listener bound in _setupViewModeSelection
    },

    // --- Core Update Functions ---

    switchSubject: function(subjectId) {
        if (!EEG_DATA.subjects || !EEG_DATA.subjects[subjectId]) {
             console.error("Subject ID not found in EEG_DATA.subjects:", subjectId); return;
        }
        console.log("Switching to subject:", subjectId);
        this.state.currentSubjectId = subjectId;
        const subjectData = EEG_DATA.subjects[subjectId];
        this.state.currentSubjectData = subjectData; // Store ref {times: [], data: [[]]}

        // Update channel count and region info from main EEG_DATA (should be static per trial)
        this.state.nAllChannels = EEG_DATA.channels?.names?.length || 0;
        this.state.allChannelRegions = EEG_DATA.channels?.ba_regions || []; // Get regions array
        if (this.state.nAllChannels === 0) { console.error("Could not determine number of channels from EEG_DATA."); }
        if (this.state.allChannelRegions.length !== this.state.nAllChannels) { console.error(`Region data length (${this.state.allChannelRegions.length}) mismatch with channel count (${this.state.nAllChannels}).`); }

        // Calculate subject-specific state
        if (!subjectData || !subjectData.times || !subjectData.data) {
             console.error("Invalid data structure for subject:", subjectId);
             this.state.nTimes = 0;
             this.state.currentSubjectMin = this.config.DEFAULT_SUBJECT_RANGE[0];
             this.state.currentSubjectMax = this.config.DEFAULT_SUBJECT_RANGE[1];
        } else {
            this.state.nTimes = subjectData.times.length;
            // Calculate overall Min/Max signal value for consistent scaling
            let overallMin = Infinity; let overallMax = -Infinity;
            subjectData.data.forEach(timeSlice => { // Iterate through time points
                 if (Array.isArray(timeSlice) && timeSlice.length === this.state.nAllChannels) {
                    timeSlice.forEach(value => { // Iterate through channels at this time point
                         if (value !== null && !isNaN(value) && isFinite(value)) {
                             if (value < overallMin) overallMin = value;
                             if (value > overallMax) overallMax = value;
                         }
                    });
                 } else { /* Optionally warn about invalid time slice format */ }
            });
            // Set min/max, handle edge cases (no data or flat data)
            if (overallMin === Infinity) { // No valid data found
                this.state.currentSubjectMin = this.config.DEFAULT_SUBJECT_RANGE[0];
                this.state.currentSubjectMax = this.config.DEFAULT_SUBJECT_RANGE[1];
            } else if (overallMin === overallMax) { // Data is flat
                const epsilon = Math.abs(overallMin * 0.1) || 0.1; // Add/subtract small value
                this.state.currentSubjectMin = overallMin - epsilon;
                this.state.currentSubjectMax = overallMax + epsilon;
            } else {
                this.state.currentSubjectMin = overallMin;
                this.state.currentSubjectMax = overallMax;
            }
        }
        console.log(`Subject ${subjectId} range: [${this.state.currentSubjectMin.toExponential(2)}, ${this.state.currentSubjectMax.toExponential(2)}] across ${this.state.nAllChannels} channels, ${this.state.nTimes} time points.`);

        // Update UI elements
        this.elements.timeSlider.max = this.state.nTimes > 0 ? this.state.nTimes - 1 : 0;
        this.elements.timeSlider.value = 0; // Reset time slider to beginning

        // Trigger subsequent updates
        this._updateActiveChannelIndices(); // Recalculate active indices based on subgroup selection
        this._setupLegend();                // Update legend based on regions potentially found in data
        this.updateTimeseriesPlot();       // Update timeseries plot (or clear if no channel selected)
        this.updateTopography();           // Redraw topography based on new subject and view mode
    },

    _updateActiveChannelIndices: function() {
        const selectedSubgroup = this.elements.subgroupSelect.value;
        // console.log("Updating active channels for subgroup:", selectedSubgroup);
        if (selectedSubgroup === "All Channels" || !EEG_DATA.channelSubgroups) {
            this.state.currentActiveChannelIndices = Array.from({ length: this.state.nAllChannels }, (_, i) => i);
        } else if (EEG_DATA.channelSubgroups[selectedSubgroup]) {
            this.state.currentActiveChannelIndices = EEG_DATA.channelSubgroups[selectedSubgroup];
        } else {
            console.warn(`Subgroup '${selectedSubgroup}' not found in EEG_DATA. Defaulting to all channels.`);
            this.state.currentActiveChannelIndices = Array.from({ length: this.state.nAllChannels }, (_, i) => i);
        }
        // console.log(`Active channel indices set: ${this.state.currentActiveChannelIndices.length}/${this.state.nAllChannels}`);
    },

    updateTopography: function() {
        if (this.state.centerX === undefined || !this.state.currentSubjectData || !this.state.channelMidpointSvgPos || this.state.channelMidpointSvgPos.length === 0) {
            console.warn("updateTopography called before initialization or no subject data/positions.");
            return;
        }
        if (!this.state.allChannelRegions || this.state.allChannelRegions.length !== this.state.nAllChannels) {
             console.error("Anatomical region data (ba_regions) is missing or length mismatch.");
             if (this.state.currentVisualizationMode === 'anatomical') {
                console.warn("Cannot display anatomical view due to missing region data.");
             }
        }

        this._updateControlDisplays();
        this._updateChannelStyles();

        // --- Mode-Specific Drawing ---
        const isHeatmapMode = this.state.currentVisualizationMode === 'heatmap';
        this.elements.heatmapGroup.style("display", isHeatmapMode ? null : "none");
        this.elements.anatomicalRegionsGroup.style("display", isHeatmapMode ? "none" : null);
        d3.select(".colorbar").style("visibility", isHeatmapMode ? "visible" : "hidden");
        d3.select("#legend-container").style("display", isHeatmapMode ? "none" : null);

        if (isHeatmapMode) {
            this._drawHeatmap();
        } else {
            this._drawAnatomicalMarkers();
        }

        // --- Common Updates (Bring interactive elements to front) ---
        this.elements.channelLinesGroup?.raise();
        this.elements.optodeDotsGroup?.raise();
        this.elements.fiducialsGroup?.raise();
        this.elements.channelClickTargetsGroup?.raise();

        this._updateCurrentTimeMarker(); // Update time marker on timeseries plot
    },

     _updateControlDisplays: function() {
        const timeIdx = parseInt(this.elements.timeSlider.value);
        const smoothingFactor = parseFloat(this.elements.smoothingSlider.value);

        // Update time display
        if (this.state.currentSubjectData?.times && timeIdx >= 0 && timeIdx < this.state.currentSubjectData.times.length) {
             this.elements.timeDisplay.textContent = this.state.currentSubjectData.times[timeIdx].toFixed(1);
        } else { this.elements.timeDisplay.textContent = timeIdx; } // Fallback to index if times array is missing
        // Update smoothing display
        this.elements.smoothingDisplay.textContent = smoothingFactor.toFixed(2);
    },

    _updateChannelStyles: function() {
        this.elements.channelLinesGroup?.selectAll('.channel-line')
            .classed('inactive', (d, i, nodes) => {
                // Get the bound datum (channel index)
                const channelIndex = d3.select(nodes[i]).datum();
                return !this.state.currentActiveChannelIndices.includes(channelIndex);
            });
            // Reset stroke color explicitly - maybe color later based on mode if desired
            // .style('stroke', '#555');
    },

    _drawHeatmap: function() {
        const timeIdx = parseInt(this.elements.timeSlider.value);
        const smoothingFactor = parseFloat(this.elements.smoothingSlider.value);

        // Check time index validity and data availability
        if (this.state.nTimes === 0 || timeIdx < 0 || !this.state.currentSubjectData?.data || timeIdx >= this.state.currentSubjectData.data.length) {
            console.warn("Invalid time index or no time data for heatmap:", timeIdx);
            this.elements.heatmapGroup.selectAll('.heat-rect').remove(); return;
        }
        const allTimeData = this.state.currentSubjectData.data[timeIdx]; // Data for all channels at this time
        if (!allTimeData || allTimeData.length !== this.state.nAllChannels) {
            console.error(`Time data length mismatch for heatmap. Expected ${this.state.nAllChannels}, got ${allTimeData?.length}`);
            this.elements.heatmapGroup.selectAll('.heat-rect').remove(); return;
        }

        // Filter signal values and 3D positions for ACTIVE channels only
        const activeValues = [];
        const activePositions3D = []; // Use DIGITIZER positions for interpolation
        this.state.currentActiveChannelIndices.forEach(index => {
            if (index < allTimeData.length && index < EEG_DATA.channels.positions3d.length) {
                activeValues.push(allTimeData[index]);
                activePositions3D.push(EEG_DATA.channels.positions3d[index]); // Get corresponding 3D pos
            } else { console.warn(`Index ${index} out of bounds during heatmap data filtering.`); }
        });

        // Update signal color scale domain and colorbar labels
        this.scales.colorScaleSignal.domain([this.state.currentSubjectMin, this.state.currentSubjectMax]);
        document.getElementById('min-val').textContent = this.state.currentSubjectMin.toExponential(1);
        document.getElementById('max-val').textContent = this.state.currentSubjectMax.toExponential(1);

        // Interpolate and Draw Heatmap Pixels
        this.elements.heatmapGroup.selectAll('.heat-rect').remove(); // Clear previous heatmap
        const pixels = [];
        const GRID_SIZE = this.config.GRID_SIZE; // Use config
        for (let i = 0; i <= GRID_SIZE; i++) {
            const svg_x = this.state.xMinGrid + i * this.state.gridStepX;
            for (let j = 0; j <= GRID_SIZE; j++) {
                const svg_y = this.state.yMinGrid + j * this.state.gridStepY;
                const targetPoint3D = this.state.gridPoints3D[i]?.[j]; // Precomputed Digitizer 3D grid point
                if (!targetPoint3D || targetPoint3D.some(isNaN)) continue; // Skip if grid point is invalid

                const value = this._interpolateValue3D(targetPoint3D, activeValues, activePositions3D, smoothingFactor);
                pixels.push({ x: svg_x, y: svg_y, value: value });
            }
        }

        this.elements.heatmapGroup.selectAll('.heat-rect')
            .data(pixels).enter().append('rect').attr('class', 'heat-rect')
            .attr('x', d => d.x - this.state.gridStepX / 2).attr('y', d => d.y - this.state.gridStepY / 2)
            .attr('width', this.state.gridStepX).attr('height', this.state.gridStepY)
            .style('fill', d => {
                if (isNaN(d.value)) return '#dddddd'; // Gray for NaN background
                 // Clamp value to the subject's min/max for consistent coloring
                const clampedValue = Math.max(this.state.currentSubjectMin, Math.min(this.state.currentSubjectMax, d.value));
                return this.scales.colorScaleSignal(clampedValue);
            });
    },

    _drawAnatomicalMarkers: function() {
        if (!this.state.allChannelRegions || this.state.allChannelRegions.length !== this.state.nAllChannels) {
             console.warn("Cannot draw anatomical regions - region data missing/invalid.");
             this.elements.anatomicalRegionsGroup.selectAll('.region-marker').remove(); return;
        }

        // Use data join based on channelMidpointSvgPos which contains {x, y, channelIndex, channelName}
        const regionMarkers = this.elements.anatomicalRegionsGroup.selectAll('.region-marker')
            .data(this.state.channelMidpointSvgPos, d => d.channelIndex); // Join by channel index

        regionMarkers.exit().remove(); // Remove markers for channels no longer present (shouldn't happen here)

        regionMarkers.enter()
            .append('circle')
            .attr('class', 'region-marker')
            .attr('r', 6) // Marker size
            .merge(regionMarkers) // Apply updates to both entering and updating elements
            .attr('cx', d => d.x) // Set/update position
            .attr('cy', d => d.y)
            .style('fill', d => {
                const regionName = this.state.allChannelRegions[d.channelIndex] || this.config.UNKNOWN_REGION_NAME;
                return this.scales.regionColorScale(regionName); // Use region color scale
             })
            .classed('inactive', d => !this.state.currentActiveChannelIndices.includes(d.channelIndex)); // Set active/inactive class
    },


    updateTimeseriesPlot: function() {
        const { selectedChannelIndex, selectedChannelName } = this.state;
        const { timeseriesSvg, tsLinePath, tsCurrentTimeMarker, tsXAxisGroup, tsYAxisGroup, timeseriesLabel } = this.elements;
        const { tsXScale, tsYScale } = this.scales;

        // Clear plot if no channel selected or no data
        if (!timeseriesSvg || !this.state.currentSubjectData || selectedChannelIndex === null) {
            tsLinePath?.attr('d', null); // Use optional chaining
            tsCurrentTimeMarker?.style('visibility', 'hidden');
            tsXAxisGroup?.selectAll('*').remove();
            tsYAxisGroup?.selectAll('*').remove();
            if (timeseriesLabel) timeseriesLabel.textContent = 'Click a channel midpoint to view its signal';
            return;
        }

        const times = this.state.currentSubjectData.times;
        let signal = [];
        try {
            // Extract signal for the selected channel index across all time points
            signal = this.state.currentSubjectData.data.map(timeSlice => timeSlice[selectedChannelIndex]);
        } catch (e) {
             console.error(`Error extracting signal for channel index ${selectedChannelIndex}`, e);
             tsLinePath?.attr('d', null);
             tsCurrentTimeMarker?.style('visibility', 'hidden');
             if (timeseriesLabel) timeseriesLabel.textContent = `Error loading data for ${selectedChannelName}`;
             return;
        }

        if (!times || times.length === 0 || !signal || signal.length !== times.length) {
             console.warn(`Inconsistent times/signal data for channel ${selectedChannelName}`);
             tsLinePath?.attr('d', null); return;
        }

        // Determine Y scale based on overall subject min/max
        let yMin = this.state.currentSubjectMin; let yMax = this.state.currentSubjectMax;
        if (yMin === yMax) { // Handle flat data
             const epsilon = Math.abs(yMin * 0.1) || 0.1; yMin -= epsilon; yMax += epsilon;
        }

        // Set scales domains
        tsXScale.domain(d3.extent(times));
        tsYScale.domain([yMin, yMax]).nice(); // Use nice() for better axis ticks

        // Draw axes
        tsXAxisGroup.call(d3.axisBottom(tsXScale).ticks(5));
        tsYAxisGroup.call(d3.axisLeft(tsYScale).ticks(5).tickFormat(d3.format(".1e"))); // Exponential format

        // Define line generator
        const line = d3.line().x(d => tsXScale(d.time)).y(d => tsYScale(d.value));
        const plotData = times.map((t, i) => ({ time: t, value: signal[i] }));

        // Draw line with transition
        tsLinePath.datum(plotData).transition().duration(50).attr("d", line);

        // Update label and time marker
        const region = this.state.allChannelRegions[selectedChannelIndex] || this.config.UNKNOWN_REGION_NAME;
        timeseriesLabel.textContent = `Signal: ${selectedChannelName} (${region})`;
        this._updateCurrentTimeMarker();
    },

     _updateCurrentTimeMarker: function() {
         const { tsCurrentTimeMarker } = this.elements;
         const { tsXScale } = this.scales;
         const { currentSubjectData, selectedChannelIndex, tsHeight } = this.state;

         if (!currentSubjectData || selectedChannelIndex === null || !currentSubjectData.times || !tsCurrentTimeMarker || !tsXScale) {
             tsCurrentTimeMarker?.style('visibility', 'hidden'); return;
         }
         const timeIdx = parseInt(this.elements.timeSlider.value);
         if (timeIdx >= 0 && timeIdx < currentSubjectData.times.length) {
             const currentTime = currentSubjectData.times[timeIdx];
             // Check if scale domain is valid before calculating position
             if (tsXScale.domain()[0] !== undefined && tsXScale.domain()[1] !== undefined) {
                  const xPos = tsXScale(currentTime);
                  // Ensure xPos is within the plot range visually
                  if (xPos >= tsXScale.range()[0] && xPos <= tsXScale.range()[1]) {
                      tsCurrentTimeMarker.attr('x1', xPos).attr('x2', xPos)
                                         .attr('y1', 0).attr('y2', tsHeight) // Use calculated height
                                         .style('visibility', 'visible');
                      return; // Exit if successful
                  }
             }
         }
         // Hide marker if time index is invalid or xPos is outside range
         tsCurrentTimeMarker.style('visibility', 'hidden');
    },

    _setupLegend: function() {
        if (!this.elements.legendContainer) return;
        this.elements.legendContainer.innerHTML = ''; // Clear previous

        // Get unique regions present in the current data (excluding None/null)
        const activeRegions = new Set(this.state.allChannelRegions.filter(r => r && r !== this.config.OUTSIDE_BRAIN_NAME));
        if (activeRegions.size === 0 && this.state.currentVisualizationMode === 'anatomical') {
            console.warn("No anatomical regions found in data to display in legend.");
            // Optionally display a message in the legend area
        }

        // Define a preferred order for the legend
        const preferredOrder = Object.keys(this.REGION_COLORS); // Use keys from defined colors

        // Filter and sort regions based on preferred order, putting others at the end
        const sortedRegions = preferredOrder.filter(region => activeRegions.has(region));
        activeRegions.forEach(region => {
            if (!preferredOrder.includes(region)) {
                sortedRegions.push(region); // Add regions found in data but not in our defined list
            }
        });

        // Setup the color scale domain and range based ONLY on the regions to be displayed
        this.scales.regionColorScale
            .domain(sortedRegions)
            .range(sortedRegions.map(r => this.REGION_COLORS[r] || this.REGION_COLORS[this.config.UNKNOWN_REGION_NAME])); // Fallback color

        // Create legend items using D3 for consistency
        const legend = d3.select(this.elements.legendContainer);
        sortedRegions.forEach(regionName => {
             const item = legend.append('div').attr('class', 'legend-item');
             item.append('div')
                 .attr('class', 'legend-color-box')
                 .style('background-color', this.scales.regionColorScale(regionName)); // Use the scale
             item.append('span').text(regionName);
        });

        // Show/hide the entire legend container based on view mode
        d3.select("#legend-container").style("display", this.state.currentVisualizationMode === 'anatomical' ? null : 'none');
    },

    // --- Event Handlers ---
    handleChannelClick: function(event, d) { // d is the bound datum {x, y, channelIndex, channelName}
         this.state.selectedChannelIndex = d.channelIndex;
         this.state.selectedChannelName = d.channelName;
         const region = this.state.allChannelRegions[d.channelIndex] || this.config.UNKNOWN_REGION_NAME;
         console.log(`Channel clicked: ${d.channelName} (Index: ${d.channelIndex}), Region: ${region}`);

         // Highlight selected channel target
         this.elements.channelClickTargetsGroup.selectAll('.channel-click-target')
             .classed('selected', datum => datum.channelIndex === this.state.selectedChannelIndex);

         this.updateTimeseriesPlot(); // Update the timeseries plot
    },


    // --- Helper Functions ---
    _interpolateValue3D: function(targetPoint3D, activeValues, activeChannelPoints3D, smoothingFactor) {
        if (!targetPoint3D || targetPoint3D.some(isNaN) || !activeChannelPoints3D || activeChannelPoints3D.length === 0) {
            return NaN;
        }
        let weightedSum = 0;
        let weightSum = 0;
        const epsilon = 1e-9;
        const power = Math.max(0.1, 2 / smoothingFactor); // Avoid zero power

        for (let i = 0; i < activeChannelPoints3D.length; i++) {
            const p3d = activeChannelPoints3D[i]; // Digitizer coords of active channel i
            const val = activeValues[i];          // Signal value of active channel i

            if (p3d === null || p3d.some(isNaN) || val === null || isNaN(val) || !isFinite(val)) {
                continue; // Skip if position or value is invalid/NaN
            }

            const dx = targetPoint3D[0] - p3d[0];
            const dy = targetPoint3D[1] - p3d[1];
            const dz = targetPoint3D[2] - p3d[2];
            const distSq = dx*dx + dy*dy + dz*dz;

            if (distSq < epsilon) { // Exact match or very close
                return val;
            }

            const distance = Math.sqrt(distSq);
            const weight = 1 / Math.pow(distance, power);

            if (isFinite(weight)) { // Avoid issues with massive weights if distance is tiny but > epsilon
                 weightedSum += val * weight;
                 weightSum += weight;
            }
        }

        if (weightSum < this.config.MIN_WEIGHT_SUM_THRESHOLD) { // Check threshold
            return NaN; // Not enough contribution from nearby channels
        }

        return weightSum > epsilon ? weightedSum / weightSum : 0; // Avoid division by zero if weightSum is tiny
    },

    _projectPointAzimuthalEquidistantJS: function(point3d_arr, center3d_arr) {
        // Ensure inputs are arrays
        if (!Array.isArray(point3d_arr) || !Array.isArray(center3d_arr) || point3d_arr.length !== 3 || center3d_arr.length !== 3) {
            console.warn("Invalid input to _projectPointAzimuthalEquidistantJS");
            return [0.0, 0.0];
        }
        const vec_x = point3d_arr[0] - center3d_arr[0];
        const vec_y = point3d_arr[1] - center3d_arr[1];
        const vec_z = point3d_arr[2] - center3d_arr[2];
        const dist = Math.sqrt(vec_x*vec_x + vec_y*vec_y + vec_z*vec_z);
        if (dist < 1e-9) return [0.0, 0.0];
        let acos_arg = vec_z / dist;
        acos_arg = Math.max(-1.0, Math.min(1.0, acos_arg)); // Clamp argument for robustness
        const theta = Math.acos(acos_arg);
        const phi = Math.atan2(vec_y, vec_x);
        const rho = theta;
        const proj_x = rho * Math.cos(phi);
        const proj_y = rho * Math.sin(phi);
        return [proj_x, proj_y];
    },

    _formatGroupName: function(name) {
        if (!name) return "";
        // Replace underscore with space, add space before uppercase letters (simple approach)
        let formatted = name.replace(/_/g, ' ');
        // Basic capitalization (first letter of each word)
        formatted = formatted.split(' ').map(word => word.charAt(0).toUpperCase() + word.slice(1).toLowerCase()).join(' ');
        return formatted.replace(/\s+/g, ' ').trim(); // Tidy spaces
    }

};

// --- Run Initialization on Page Load ---
document.addEventListener('DOMContentLoaded', () => {
    // Check if EEG_DATA was loaded (basic check, more detailed in init)
    if (typeof EEG_DATA !== 'undefined' && Object.keys(EEG_DATA).length > 0) {
        fnirsVisApp.init();
    } else {
        // Handle the case where EEG_DATA is missing immediately
        console.error("DOMContentLoaded: EEG_DATA is missing or empty. Cannot initialize visualization.");
        fnirsVisApp._showFatalError("Visualization data (EEG_DATA) was not loaded correctly.");
    }
});