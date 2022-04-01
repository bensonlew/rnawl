/* globals $: true, jQuery: true, d3: true,  window: true, document: true */
/* Author: ChengXihui */

; (function($, window, document) {
    "use strict";

    $.extend({
        // fnName: 函数名; el: 容器id; contents: 传入的数据
        d3chart: function(fnName, el, contents) {
            var fn = $.extend({}, D3charts), handler;
            if (fn[fnName] === undefined) {
                console.log("没有这个函数： " + fnName);
                return false;
            };
            handler = fn[fnName](el, contents);
            return handler;
        }
    });


    var D3charts = {

        nPie: function (ele, contents){
                var opts = {
                    data: [],
                    title: "",
                    margin: {
                        top: 80,
                        right: 75,
                        bottom: 40,
                        left: 40
                    },

                    // padding between pie chart groups
                    outerPadding: 120,

                    // padding between pie and its legend
                    innerPadding: 40,

                    // donutRatio = outerRadius / innerRadius
                    donutRatio: 0.35,
                    groupName: ["BP", "CC", "MF"],
                    size: { width: 1500, height: 200 }
                };

                // overwrite default opts
                Object.keys(contents).forEach(function(key) {
                    if (key in opts) {
                        opts[key] = contents[key];
                    }
                });
                var svgWidth = opts.size.width,
                    svgHeight = opts.size.height,
                    margin = opts.margin,

                    outerPadding = opts.outerPadding,
                    innerPadding = opts.innerPadding,

                    plotWidth = opts.size.width - opts.margin.left - opts.margin.right,

                    // legend rect box size = legendRectSize X legendRectSize
                    legendRectSize = 12,

                    // pie chart radius
                    outerRadius,

                    isVisible = true;


                var colors = d3.scale.category20();

                var svg = d3.select("#" + ele).html('').append('svg').attr({
                    width         : svgWidth,
                    height        : svgHeight,
                    version       : "1.1",
                    xmlns         : "http://www.w3.org/2000/svg",
                    "font-size"   : 12,
                    "font-family" : "Arial"
                });

                // testing text size at the fixed font size and font face
                function fetchTextWidth(svg, text){
                    var testNode = svg.append("text").text(text).node();
                    var f;
                    try{
                        f = testNode.getBBox().width;
                        if(f === 0){ f = text.lenght * 8; isVisible = false;}
                    }
                    catch(e){
                        f = text.lenght * 8;
                    }finally{
                        testNode.remove();
                    }
                    return f;
                }

                // data processing
                var myData = opts.data.map(function(pie){
                    var values = [], names = [], metas = [], legendMaxWidth = 0, x = 0;
                    var nameTextLenght;
                    pie.forEach(function(obj){
                        values.push(obj.value);
                        names.push(obj.name);
                        metas.push(obj.meta);
                    });

                    names.forEach(function(name){
                        nameTextLenght = fetchTextWidth(svg, name);
                        if(legendMaxWidth < nameTextLenght){ legendMaxWidth = nameTextLenght; }
                    });
                    return {
                        "values"          : values,
                        "names"           : names,
                        "metas"           : metas,
                        "legendMaxWidth"  : legendMaxWidth + legendRectSize + 2,
                        "x"               : x
                    };
                });


                function calculateRadius(){
                    var legends = myData.map(function(d){ return d.legendMaxWidth; });
                    var legendTotalWidth = legends.reduce(function(sum, value) { return sum + value; }, 0);
                    var diameter = (plotWidth - (innerPadding * (myData.length) + outerPadding * (myData.length - 1)) - legendTotalWidth) / 3;
                    if (diameter < 120){
                        svgWidth += (120 - diameter) * 3;
                        plotWidth = opts.size.width - opts.margin.left - opts.margin.right;
                        svg.attr("width", svgWidth);
                        diameter = 120;
                    }
                    outerRadius = diameter / 2;
                }
                calculateRadius();
                myData.forEach(function(d, i){
                    if(i === 0){
                        d.x = 0; d.deltaX = outerRadius * 2 + innerPadding + d.legendMaxWidth + outerPadding;
                    }else{
                        d.x = myData[i - 1].deltaX;
                        d.deltaX = d.x + outerRadius * 2 + innerPadding + d.legendMaxWidth + outerPadding;
                    }
                });

                function onMouseover(){
                    d3.select(this).attr({
                        stroke: '#666',
                        "stroke-width": 2
                    });
                    var index = d3.select(this).datum().index + 1;
                    var ele = d3.select(this.parentNode.parentNode.nextElementSibling.querySelector("g:nth-child(" + index + ")"));
                    ele.style({
                        "font-weight": "bold",
                        fill: colors.range()[index - 1]
                    });
                }
                function onMouseout(){
                    d3.select(this).attr({
                        stroke: 'none',
                        "stroke-width": 0
                    });
                    var index = d3.select(this).datum().index + 1;
                    var ele = d3.select(this.parentNode.parentNode.nextElementSibling.querySelector("g:nth-child(" + index + ")"));
                    ele.style({
                        "font-weight": "normal",
                        fill: "#000"
                    });
                }

                var mainPlot = svg.append('g')
                    .attr('class', 'pie-main-plot')
                    .attr('transform', 'translate(' + [margin.left, margin.top] + ')');

                var pieGroup = mainPlot.selectAll('g.single-pie')
                    .data(myData)
                    .enter()
                    .append('g')
                    .attr({
                        class: 'single-group-pie',
                        transform: function(d){
                            return 'translate(' + [d.x, 0] + ')';
                        }
                    });
                var arcGenerator = d3.svg.arc()
                    .outerRadius(outerRadius)
                    .innerRadius(outerRadius * opts.donutRatio)
                    .startAngle(function(d) {return d.startAngle;})
                    .endAngle(function(d) {return d.endAngle;});

                var piePlot = pieGroup.append('g')
                    .attr('class', 'pie-plot')
                    .attr("transform", "translate(" + [outerRadius, outerRadius + 50] + ")")
                    .selectAll('g')
                    .data(function(d){ return d3.layout.pie().sort(null)(d.values); })
                    .enter()
                    .append("g");
                piePlot.append('path').attr({
                    class: 'pie-segement',
                    d: arcGenerator,
                    fill: function(d, i) {
                        return colors.range()[i];
                    },
                }).datum(function(d, i) {
                    var meta = d3.select(this.parentNode.parentNode).datum().metas[i];
                    d.meta = meta;
                    d.index = i;
                    return d;
                })
                .on('mouseover', onMouseover)
                .on('mouseout', onMouseout);

                piePlot.append('g')
                    .attr('transform',function(d) {
                    var angle = (d.startAngle + d.endAngle) / 2 * 180 / Math.PI - 90;
                    return "rotate(" + angle + ")" + 'translate(' + [outerRadius, 0] + ')';
                })
                .append('text').text(function(d) {
                    return d.value;
                })
                .attr({
                    x: 8,
                    dy: '0.55em',
                    transform: function (d) {
                        var textAngle = d.startAngle + (d.endAngle - d.startAngle) * 0.5;
                        return textAngle > Math.PI ? 'rotate(180)translate(-16)' : null;
                    },
                    'text-anchor': function (d) {
                        var textAngle = d.startAngle + (d.endAngle - d.startAngle) * 0.5;
                        return textAngle > Math.PI ? 'end' : null;
                    }
                });


                var pieLegend = pieGroup.append('g')
                    .attr('class', 'legend pie-legend')
                    .attr('transform', 'translate(' + [2 * outerRadius + innerPadding, 0] + ')');

                var legendItem = pieLegend.selectAll('g')
                    .data(function(d){ return d.names; })
                    .enter()
                    .append('g')
                    .attr('transform', function(d,i){return 'translate(' + [0, i * 20] + ')';});
                legendItem
                    .append('rect')
                    .attr({
                        width         : legendRectSize,
                        height        : legendRectSize,
                        fill          : function(d,i){ return colors.range()[i]; },
                        stroke        : 'none',
                        "stroke-width": 0
                    });
                legendItem
                    .append('text')
                    .text(function(d){ return d; })
                    .attr({
                        transform: 'translate(' + [17, 10] + ')',
                    });

                // add group name text at the top left corner of every single chart
                svg.selectAll('.single-group-pie')
                .each(function(d, i){
                    d3.select(this)
                        .append('text')
                        .text(opts.groupName[i])
                        .attr('font-size', 15)
                        .style('font-weight', 'bold');
                });

                if (opts.title) {
                    svg.append('g').attr({
                        transform: 'translate(' + [svgWidth / 2, 25] + ')',
                        class: 'svg pie-main-title',
                    }).append('text').text(opts.title).attr({
                        "font-size": 20,
                        "text-anchor": "middle",
                    });
                }

                function noConflict(node1, node2, span){
                    if(!span || span < 10){ span = 10; }
                    var box1 = node1.getBoundingClientRect(), box2 = node2.getBoundingClientRect();
                    var trans2 = node2.getAttribute("transform");
                    trans2 = d3.transform(trans2).translate;
                    if((box2.left - box1.right) < 0){
                        node2.setAttribute("transform", "translate(" + [(box1.right - box2.left + trans2[0]) + span, trans2[1]]+ ")");
                    }
                }

                // if the graphic is visible using getBBox to calculate max height and reSize svg height attribute
                var maxGroupHeight = 0;
                if(isVisible){
                    svg.selectAll(".single-group-pie").each(function(){
                        var temp;
                        temp = this.getBBox().height;
                        if(maxGroupHeight < temp){ maxGroupHeight = temp; }
                    });
                    svg.attr("height", maxGroupHeight + margin.top + margin.bottom);

                    pieGroup.each(function(){
                        noConflict(this.querySelector(".pie-plot"), this.querySelector(".pie-legend"), innerPadding);
                    });

                    pieGroup[0].forEach(function(ele){
                        if(ele.nextElementSibling){
                            noConflict(ele, ele.nextElementSibling, outerPadding);
                        }
                    });

                    var mainPlotRight = svg.select(".pie-main-plot").node().getBoundingClientRect().right,
                        svgBoxRight = svg.node().getBoundingClientRect().right;
                    if(mainPlotRight > svgBoxRight){
                        svg.attr("width", svgWidth + (mainPlotRight - svgBoxRight) + 20);
                    }
                }

        },
        violin: function (ele, content){
            var o = {
                title       : "",
                margin      : {top: 50, right: 50, bottom: 50, left: 55},
                paddingRatio: 0.2,
                size        : {width: 600, height: 400},
                bins        : 20,
                show_legend : false,
                show_plotbox: false,
                x_label     : '',
                y_label     : '',

                // linear, step, step-before, step-after, basis,
                // basis-open, cardinal, cardinal-open, monotone
                // only "step-after" is appropriate from the statistics perspective.
                interpolate : 'step-after',
                data        : []
            };

            Object.keys(content).map(function(e){
                if(e in o){
                    o[e] = content[e];
                }else{
                     console.log('d3violin: Invalid option name:\t' + e + '\n');
                }
            });

            if(!['step-after', 'basis'].includes(o.interpolate)){
                throw new Error('d3violin: interpolate type is WRONG!');
            }
            if(!o.data.length){
                throw new Error('d3violin: Data is EMPTY!');
            }

            // data processing
            var groupQty = o.data.length,
                maxValue = Math.ceil( d3.max(o.data, function(e, i){ return d3.max(e.data); }) ),
                minValue = Math.floor( d3.min(o.data, function(e, i){ return d3.min(e.data); }) );
            var tickNamesLength = 0;
            var tickNames = o.data.map(function(e){
                tickNamesLength = Math.max(e.name.length, tickNamesLength);
                return e.name;
            });

            if(new Set(tickNames).size < tickNames.length){
                throw new Error('d3violin: Same Name OCCURRENCES: ' + tickNames);
            }

            if(tickNamesLength > 7 || tickNames.length > 12) {
                o.margin.bottom += 20;
                // o.margin.left += 10;
            }

            var plotWidth  = o.size.width - o.margin.left - o.margin.right,
                plotHeight = o.size.height - o.margin.top - o.margin.bottom;

            var tooltip = d3.select("body").append("div")
                .attr("class","violin_tooltip") //用于css设置类样式
                .attr("opacity",0.0);


            // define scale and axis
            var xScale = d3.scale.ordinal()
                    .domain(tickNames)
                    .rangeBands([0, plotWidth], o.paddingRatio, o.paddingRatio),

                yScale = d3.scale.linear()
                    .domain([minValue, maxValue])
                    .range([0, plotHeight]),

                xAxis = d3.svg.axis()
                    .orient('bottom')
                    .outerTickSize(0)
                    .scale(xScale),

                yAxis = d3.svg.axis()
                    .orient('left')
                    .scale(yScale.copy().domain([maxValue, minValue]));

            // create svg and divided into three parts
            var svg = d3.select("#" + ele)
                .html('')
                .append('svg')
                .attr({
                    width         : o.size.width,
                    height        : o.size.height,
                    version       : "1.1",
                    xmlns         : "http://www.w3.org/2000/svg",
                    class         : "svg-violin",
                    "font-size"   : "12px",
                    "font-family" : "Arial",
                });

                var plotMain = svg.append('g')
                    .attr('class', 'plot-main')
                    .attr('transform', 'translate(' +
                          [o.margin.left, o.margin.top] + ')');

            svg.append('g')
                .attr({
                    transform: 'translate(' +[o.margin.left, o.margin.top] + ')',
                    class: 'axis y-axis'
                })
                .call(yAxis)
                .selectAll('path, line')
                .attr({
                    fill          : 'none',
                    stroke        : '#000',
                    'stroke-width': 0.5
                });
            svg.append('g')
                .attr({
                    transform: 'translate(' +
                        [o.margin.left, o.margin.top + plotHeight] + ')',
                    class: 'axis x-axis'
                })
                .call(xAxis)
                .call(function(e){

                    e.selectAll("line")
                        .attr({
                            fill: 'none',
                            stroke: 'transparent'
                        });
                    e.select("path.domain")
                        .attr({
                            fill: 'none',
                            stroke: "#000",
                            "stroke-width": 0.5
                        });

                    if(tickNamesLength > 7) {
                        e.selectAll("text").each(function(d) {
                            var ele = d3.select(this);
                            var l = this.getComputedTextLength();
                            ele.attr("text-anchor", "start")
                            .attr("transform", `translate(-${l/4},${l/2}),rotate(-60)`);
                        });
                    } else if(tickNames.length > 12) {
                        e.selectAll("text").each(function(d) {
                            var ele = d3.select(this);
                            var l = this.getComputedTextLength();
                            ele.attr("text-anchor", "start")
                            .attr("transform", `translate(-${l/2},${l/2}),rotate(-60)`);
                        });
                    } else {

                    }

                });



            // prepare for drawing
            var boxWidth = xScale.rangeBand();
            var violinGroups = plotMain.selectAll('g.violin-group')
                .data(o.data)
                .enter()
                .append('g')
                .each(function(d,i){
                    d3.select(this).attr({
                        class    : "violin-group",
                        transform: 'translate(' + [
                            xScale.range()[i] + 0.5 * boxWidth,
                            plotHeight
                        ] + ')',
                    });
                });

            function addViolin(selection, hData, hArea, hline){
                var halfViolin = selection.append('g')
                    .attr('fill', function(d){return d.color;});

                halfViolin.append('path')
                    .datum(hData)
                    .attr('class', 'path-area')
                    .attr('d', hArea);

                halfViolin.append('path')
                    .datum(hData)
                    .attr('d', hline)
                    .attr('class', 'path-line')
                    .attr('fill', 'none');

                var cloneNode = halfViolin.node().cloneNode(true);
                cloneNode = d3.select(selection.node().insertBefore(cloneNode, null));
                cloneNode.attr('transform', 'rotate(-90) scale(1,-1)');
                halfViolin.attr('transform', 'rotate(-90)');
            }
            function addPlotBox(selection, hData, boxWidth) {

                // Quantile Probability
                // https://en.wikipedia.org/wiki/Quantile
                var probs = [0.05, 0.25, 0.5, 0.75, 0.95];
                probs = probs.map(function(e){
                    return yScale(d3.quantile(hData, e));
                });
                var mean = yScale(d3.mean(hData));
                var stickData = {
                        w1: probs[0],
                        q1: probs[1],
                        q2: probs[2],
                        q3: probs[3],
                        w2: probs[4],
                        mean: mean
                    };

                var plotBox = selection.append('g')
                    .datum(stickData)
                    .attr('opacity', 0.5)
                    .style('cursor', 'pointer')
                    .attr('class', 'violin-plotbox');

                var width = Math.max(10, boxWidth * 0.25);

                // add two whiskers and a center line
                plotBox.append('line')
                    .attr({
                        x1    : -0.5 * width,
                        x2    : 0.5 * width,
                        y1    : -probs[4],
                        y2    : -probs[4],
                        class : 'plot-box-line'
                    });
                plotBox.append('line')
                    .attr({
                        x1    : -0.5 * width,
                        x2    : 0.5 * width,
                        y1    : -probs[0],
                        y2    : -probs[0],
                        class : 'plot-box-line'
                    });
                plotBox.append('line')
                    .attr({
                        x1    : 0,
                        x2    : 0,
                        y1    : -probs[0],
                        y2    : -probs[4],
                        class : 'plot-box-line'
                    });

                // median line, box and mean circles
                plotBox.append('rect')
                    .attr({
                        x             : -0.5 * width,
                        y             : -probs[3],
                        width         : width,
                        height        : probs[3] - probs[1],
                        stroke        : '#000',
                        'stroke-width': 1,
                        fill          : '#000'
                    });
                plotBox.append('line')
                    .attr({
                        x1            : -0.5 * width,
                        x2            : 0.5 * width,
                        y1            : -probs[2],
                        y2            : -probs[2],
                        'stroke-width': 2,
                        stroke        : '#fff'
                    });
                plotBox.append('circle')
                    .attr({
                        cx  : 0,
                        cy  : -mean,
                        r   : Math.max(4, 0.4 * width),
                        fill: '#fff'
                    });
                plotBox.append('circle')
                    .attr({
                        cx  : 0,
                        cy  : -mean,
                        r   : Math.max(4, 0.4 * width) - 2,
                        fill: '#000'
                    });

                // line attributes settings
                plotBox.selectAll('line.plot-box-line').attr({
                    fill: 'none',
                    stroke: '#000',
                    'stroke-width': 1
                });

            }

            function mouseover(){
                d3.select(this.parentNode).attr('opacity', 1);

                var datum = d3.select(this).datum();
                var self = {};
                for(var e in datum){
                    self[e] = Math.round(datum[e] * 100) / 100;
                };


                var page_x     = d3.event.pageX;
                var page_y     = d3.event.pageY+20;

                var ul = '<ul style="margin:0;padding:3px;list-style:none;border:1px dashed #555">' +
                    '<li style="margin:0;padding:0;">' + '<span style="font-weight:bolder;display:inline-block;width:35px;">W2: </span>' + self.w2 + '</li>' +
                    '<li style="margin:0;padding:0;">' + '<span style="font-weight:bolder;display:inline-block;width:35px;">Q3: </span>' + self.q3 + '</li>' +
                    '<li style="margin:0;padding:0;">' + '<span style="font-weight:bolder;display:inline-block;width:35px;">Q2: </span>' + self.q2 + '</li>' +
                    '<li style="margin:0;padding:0;">' + '<span style="font-weight:bolder;display:inline-block;width:35px;">Q1: </span>' + self.q1 + '</li>' +
                    '<li style="margin:0;padding:0;">' + '<span style="font-weight:bolder;display:inline-block;width:35px;">W1: </span>' + self.w1 + '</li>' +
                    '<li style="margin:0;padding:0;border-top:1px dashed #555;">' + '<span style="font-weight:bolder;display:inline-block;width:35px;">Mean: </span>' + self.mean + '</li>' +
                    '</ul>';

                tooltip.html('<b>W2:'+self.w2+"</b><br/><b>Q3:"+self.q3 + "</b><br/><b>Q2:" +self.q2 + "</b><br/><b>Q1:" + self.q1
                        + "</b><br/><b>W1:" + self.w1 + "</b><br/><b>Mean:" + self.mean + "</b>")
                    .style("position", "absolute")
                    .style("left",page_x+"px")
                    .style("top",page_y+"px")
                    .style("opacity",0.9)
                    .style('padding', '5px');

                // var ul = '<ul style="margin:0;padding:3px;list-style:none;border:1px dashed #555">' +
                //     '<li style="margin:0;padding:0;">' + '<span style="font-weight:bolder;display:inline-block;width:35px;">W2: </span>' + self.w2 + '</li>' +
                //     '<li style="margin:0;padding:0;">' + '<span style="font-weight:bolder;display:inline-block;width:35px;">Q3: </span>' + self.q3 + '</li>' +
                //     '<li style="margin:0;padding:0;">' + '<span style="font-weight:bolder;display:inline-block;width:35px;">Q2: </span>' + self.q2 + '</li>' +
                //     '<li style="margin:0;padding:0;">' + '<span style="font-weight:bolder;display:inline-block;width:35px;">Q1: </span>' + self.q1 + '</li>' +
                //     '<li style="margin:0;padding:0;">' + '<span style="font-weight:bolder;display:inline-block;width:35px;">W1: </span>' + self.w1 + '</li>' +
                //     '<li style="margin:0;padding:0;border-top:1px dashed #555;">' + '<span style="font-weight:bolder;display:inline-block;width:35px;">Mean: </span>' + self.mean + '</li>' +
                //     '</ul>';
                // d3.select('#' + ele)
                //     .select('.svg-tooltips')
                //     .html(ul)
                //     .style('display', 'block')
                //     .style('left', (d3.event.pageX + 5) + 'px')
                //     .style('top', (d3.event.pageY - 15) + 'px');
            }
            function mouseout(){
                tooltip.style('opacity', 0);
                // d3.select(this.parentNode).attr('opacity', 0.5);

                // d3.select('#' + ele)
                //     .select('.svg-tooltips')
                //     .style('display', 'none');
            }


            for(var i = 0; i < o.data.length; i++){
                var histogram = d3.layout.histogram()
                        .bins(o.bins);

                o.data[i].data = o.data[i].data.sort(d3.ascending);

                var histData = histogram(o.data[i].data),

                    xBoxScale = d3.scale.linear()
                        .domain([0, d3.max(histData.map(function(e){return e.y;}))])
                        .range([0, 0.5 * boxWidth]),

                    area = d3.svg.area()
                        .x(function(d){return yScale(d.x);})
                        .y(function(d){return xBoxScale(d.y);})
                        .interpolate(o.interpolate)
                        .y0(0),

                    line = d3.svg.line()
                        .x(function(d){return yScale(d.x);})
                        .y(function(d){return xBoxScale(d.y);})
                        .interpolate(o.interpolate);

                addViolin(d3.select(violinGroups[0][i]), histData, area, line);

                if(o.show_plotbox){
                    addPlotBox(d3.select(violinGroups[0][i]), o.data[i].data, boxWidth);

                    d3.select(svg.node().parentNode)
                        .style('position', 'relative')
                        .append('div')
                        .attr('class', 'svg-tooltips')
                        .style({
                            position          : 'absolute',
                            display           : "none",
                            'min-width'       : "30px",
                            margin            : '0',
                            padding           : '0',
                            'font-size'       : '10px',
                            'background-color': 'rgba(200, 200, 200, 0.5)'
                        });

                    svg.selectAll('g.violin-plotbox *')
                        .on('mouseover', mouseover)
                        .on('mouseout', mouseout);
                }
            };

            // add legend if needed
            if(o.show_legend){
                var span = 20;

                var legend = svg.append('g')
                    .attr('transform', 'translate(' +
                          [o.size.width - o.margin.right, o.margin.top] + ')')
                    .attr('class', 'legend');

                var gLegend = legend.selectAll('g')
                    .data(o.data.map(function(d,i){
                            return {
                                color : d.color,
                                name  : d.name
                            };
                        }))
                    .enter()
                    .append('g')
                    .attr('transform', function(d, i){
                        return 'translate(' + [0, i * span] + ')';
                    });

                gLegend.append('rect')
                    .attr({
                        width         : span - 5,
                        height        : span - 5,
                        fill          : function(d){return d.color;},
                        stroke        : 'none',
                        'stroke-width': 0
                    });
                gLegend.append('text')
                    .attr({
                        x: span,
                        y: span - 8,
                    })
                    .text(function(d){return d.name;});
            } // end if

            // add title if needed
            if(o.title){
                svg.append('text')
                .attr({
                    x             : 0.5 * o.size.width,
                    y             : 25,
                    'text-anchor' : 'middle',
                    class         : 'svg-title'
                })
                .style('font-size', '20px')
                .text(o.title);
            };

            // add x_label & y_label
            if(o.x_label){
                svg.append('text')
                .attr({
                    x             : o.margin.left + 0.5 * (o.size.width - o.margin.left),
                    y             : o.size.height - 10,
                    'text-anchor' : 'middle',
                    class         : 'svg-label'
                })
                .style('font-size', '15px')
                .text(o.x_label);
            }
            if(o.y_label){
                svg.append('text')
                .attr({
                    x             : 15,
                    y             : o.margin.top + 0.5 * plotHeight,
                    'text-anchor' : 'middle',
                    'transform'   : 'rotate(-90, ' + [15, o.margin.top + 0.5 * plotHeight] + ')',
                    class         : 'svg-label'
                })
                .style('font-size', '15px')
                .text(o.y_label);
            }
        },
        pieDistribute: function(ele, content){
            var opt = {
                title: 'Chrome Distribution',
                margin: {top: 90, right: 90, bottom: 90, left: 90},

                // radiusRatio: [r1, r2, r3]
                // [r1, r2] -> inner donut's innerRadius and outerRadius
                // [r2, r3] -> two dounts' padding
                // [r3] -> outer donut's innerRadius
                radiusRatio: [0.5, 0.9, 0.95],

                // if assignment, this value will take precedence
                tickProportion: undefined,

                // out donut padAngle
                padAngle: 0.03,
                data: [],
                size: {width: 500, height: 500},
            };

            Object.keys(content).map(function(e){
                if(e in opt){
                    opt[e] = content[e];
                }else{
                     console.log('pieDistribute: Invalid option name:\t' + e + '\n');
                }
            });

            if(!opt.data.length){
                throw new Error('pieDistribute: Data is EMPTY!');
            }

            // -> [[a, c, d, ...], [ ... ], [ ... ], ....]
            var temp = opt.data.map(function (D) {
                    return D.children.map(function (d) {
                        return d[1];
                    });
                }),

                // children max & min values for radiusScale
                extent = d3.extent(Array.prototype.concat.apply([], temp)),

                // tick proportion of d.value for ticks calculation
                tickProportion = d3.sum(opt.data.map(function (d) {
                        return d.value;
                    })) / 360;

            // sort children ascendantly according to [0][0]
            // the sorted value will be used for inner pie generation
            opt.data.map(function (e) {
                e.children = e.children.sort(function (a, b) {
                    return a[0][0] - b[0][0];
                });
            });

            var plotWidth  = opt.size.width - opt.margin.left - opt.margin.right,
                plotHeight = opt.size.height - opt.margin.top - opt.margin.bottom,
                maxRadius  = Math.min(plotWidth, plotHeight) * 0.5;

            var radiusScale = d3.scale.linear()
                .domain(extent)
                .range([maxRadius * opt.radiusRatio[0],
                        maxRadius * opt.radiusRatio[1]]);

            var arc = d3.svg.arc()
                .innerRadius(maxRadius * opt.radiusRatio[2])
                .outerRadius(maxRadius),

                innerArc = d3.svg.arc()
                .innerRadius(maxRadius * opt.radiusRatio[0]);

            var pie = d3.layout.pie()
                .sort(null)
                .padAngle(opt.padAngle)
                .value(function(d){return d.value;})(opt.data),

                innerPie = d3.layout.pie();

            // @ plot begin
            var svg = d3.select('#' + ele)
                .append('svg')
                .attr({
                    width   : opt.size.width,
                    height  : opt.size.height,
                    version : "1.1",
                    xmlns   : "http://www.w3.org/2000/svg"
                })
                .style({
                    "font-size"   : "10px",
                    "font-family" : "Arial",
                });
            var gCenter = svg.append('g')
                .attr({
                    class: 'gCenter',
                    transform: 'translate(' +
                        [0.5 * opt.size.width, 0.5 * opt.size.height] + ')'
                });

            var gPie = gCenter.selectAll('g.gPie')
                .data(pie)
                .enter()
                .append('g')
                .attr('class', 'gPie');

            gPie.append('path')
                .attr({
                    class: 'outer-donut',
                    fill: function(d){return d.data.color;},
                    d: arc
                });

            gPie.each(function(D, index){
                var self = d3.select(this);
                var pieData = D.data.children.map(function (d) {
                        return d[0][1] - d[0][0];
                    }),
                    pieRadius = D.data.children.map(function (d) {
                        return radiusScale(d[1]);
                    });

                innerPie
                    .sort(null) // the data is pre-sorted
                    .startAngle(D.startAngle + (0.5 * D.padAngle))
                    .endAngle(D.endAngle - (0.5 * D.padAngle));

                pieData = innerPie(pieData);
                pieRadius.map(function(d, i){
                    pieData[i].outerRadius = pieRadius[i];
                });

                // generate inner donut
                self.append('g')
                    .selectAll('.inner-donut')
                    .data(pieData)
                    .enter()
                    .append('path')
                    .attr({
                        class : 'inner-donut',
                        fill  : D.data.color,
                        d     : innerArc
                    });

                // generate ticks and labels outside the outer donut
                var gTicks = self.append('g')
                    .attr('class', 'ticks')
                    .selectAll('g')
                    .data(groupTicks)
                    .enter()
                    .append('g')
                    .attr({
                        transform: function(d){
                            return 'rotate(' + (d.angle * 180 / Math.PI - 90) +')translate(' + maxRadius +', 0)';
                        }
                    });

                gTicks.append('line')
                    .attr({
                        x1: 0, y1: 0,
                        x2: 5, y2: 0,
                        stroke        : '#000',
                        'stroke-width': function(d, i){return i % 5 ? 0.5 : 1;}
                    });

                gTicks.append('text')
                    .attr({
                        x: 8,
                        dy: '0.55em',
                        transform: function (d) {
                            return d.angle > Math.PI ?'rotate(180)translate(-16)' : null;
                        },
                        'text-anchor': function (d) {
                            return d.angle > Math.PI ? 'end' : null;
                        }
                    })
                    .text(function (d) {
                        return d.label;
                    });

                // section names
                var textAngle = D.startAngle + (D.endAngle - D.startAngle) * 0.5;
                self.append('g')
                    .attr('transform',
                          'rotate(' + (textAngle * 180 / Math.PI - 90) + ')' +
                          'translate(' + [maxRadius + 25, 0] + ')')
                    .append('text')
                    .attr({
                        transform     : textAngle > Math.PI ? 'rotate(180)' : null,
                        'text-anchor' : textAngle > Math.PI ? 'end' : null,
                        dy            : 8
                    })
                    .style('font-size', '12px')
                    .text(D.data.name);

            });

            // information center
            gCenter.append('g')
                .attr('class', 'svg-information')
                .call(function(sel){
                    // add title to center of the chart
                    if(opt.title){
                        sel.append('text')
                            .attr({
                                class: 'title',
                                'text-anchor': 'middle'
                            })
                            .style('font-size', '12px')
                            .text(opt.title);
                    }
                    var minInterval = opt.tickProportion ? opt.tickProportion : round(tickProportion);

                    // add ticks interval information
                    sel.append('text')
                        .attr({
                            class         : 'ticks-Interval',
                            'text-anchor' : 'middle',
                            y             : 15
                        })
                        .text('min unit: ' + d3.format('e')(minInterval.toString()));

                    // when hovering on a specific group,
                    // this text will show group value
                    sel.append('text')
                        .attr({
                            class         : 'gPie-value',
                            fill          : 'red',
                            'text-anchor' : 'middle',
                            y             : 30
                        })
                        .style('opacity', 0)
                        .text('...');

                    }); // end call

            // events
            gPie
                .on('mousemove', function(){
                    d3.event.stopPropagation();

                    d3.select(this.parentNode).selectAll('g.gPie').style('opacity', 0.5);
                    d3.select(this).style('opacity', 1);

                    gCenter.select('text.gPie-value')
                        .text('[ ' + d3.select(this).datum().value + ' ]')
                        .style('opacity', 1);

                    d3.select(this).select('.referance-line').remove();
                    var coordinats = d3.mouse(this),
                        rad = Math.atan2(-coordinats[1], coordinats[0]);

                    // add a reference line on hovering
                    d3.select(this).append('line')
                        .attr({
                            class             : 'referance-line',
                            stroke            : '#000',
                            'stroke-dasharray': '5,5',
                            'stroke-width'    : 1,
                            x2                : maxRadius * Math.cos(rad),
                            y2                : maxRadius * Math.sin(rad) * -1
                        });
                })
                .on('mouseout', function(){
                    d3.select(this.parentNode)
                        .selectAll('g.gPie')
                        .style('opacity', 1);

                    gCenter.select('text.gPie-value')
                        .style('opacity', 0);

                    d3.select(this).select('.referance-line').remove();
                });

            function round(nub){
                if(nub >= 1){
                    nub = Math.ceil(nub);
                    var length = nub.toString().length;
                    var e = '1e' + (length - 1);

                    return Math.round(nub / Number(e)) * Number(e);
                }else{
                    // eg: '0.000546'
                    var deep = nub.toString().search(/[1-9]/);
                    var e = '1e' + (deep - 1);

                    return Math.round(nub * Number(e)) / Number(e);
                }
            }

            function groupTicks(d) {
                var k = (d.endAngle - d.startAngle - d.padAngle) / d.value;
                var proportion = opt.tickProportion ?
                    opt.tickProportion : round(tickProportion);
                return d3.range(0, d.value, proportion).map(function (v, i) {
                    return {
                        angle: v * k + d.startAngle + 0.5 * d.padAngle,
                        label: i % 5 ? null : v / proportion
                    };
                });
            }
        },
    };

    //调用方法
    //$.d3chart("nPie",el,contents)

})(jQuery, window, document);