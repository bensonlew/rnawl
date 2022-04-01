/* globals $: true, jQuery: true, jvenn: true */

jQuery.jvenn_report = {
    showVenn:function(container, contents){
       function Colorlist_change(StartColorList){
            var colorArr = [];
            for(var i=0;i<StartColorList.length;i++){
              var newcolor = this.colorRgb(StartColorList[i]);
              colorArr.push(newcolor);
            }
            return colorArr;
      };

      Colorlist_change.prototype.colorRgb = function (sColor) {
          var reg = /^#([0-9a-fA-f]{3}|[0-9a-fA-f]{6})$/;
          sColor = sColor.toLowerCase();
          if (sColor && reg.test(sColor)) {
              // 处理简写16进制颜色: #af3 => #aaff33
              if (sColor.length === 4) {
                  var sColorNew = "#";
                  for (var i = 1; i < 4; i += 1) {
                      sColorNew += sColor.slice(i, i + 1).concat(sColor.slice(i, i + 1));
                  }
                  sColor = sColorNew;
              }
              // 处理完整的16进制颜色： #a8f6c3 => rgb(r,g,b)
              var sColorChange = [];
              for (var i = 1; i < 7; i += 2) {
                  sColorChange.push(parseInt("0x" + sColor.slice(i, i + 2)));
              }
              var sColorChangeRGB = 'rgb(' + sColorChange + ')';
              return sColorChangeRGB;
          } else {
              // 如果不属于上述两种直接返还原值，但这个估计不灵，因为jvenn要求只能rgb模式
              return sColor;
          }
      };

      var newcolorlist = new Colorlist_change(contents.params.color);
      $("#"+container).jvenn({
        series: contents.data,
        colors: newcolorlist,
        fontSize: contents.params.fontSize,
        fontFamily: contents.params.fontFamily,
        searchInput:  $("#"+ contents.params.div_search_field),
        searchStatus: $("#"+contents.params.div_search_status),
        displayMode: 'classic',
        shortNumber: false,
        displayStat1: contents.params.displayStat1,
        displayStat2: contents.params.displayStat2,

        fnClickCallback: function() {
            var value = "";
            if (this.listnames.length == 1) {
                value += "Elements only in ";
            } else {
                value += "Common elements in ";
            }
            for (var name in this.listnames) {
                value += this.listnames[name] + " ";
            }
            value += ":\n";
            for (var val in this.list) {
                value += this.list[val] + "\n";
            }
            $("#" + contents.params.div_id_names).val(value);
        }
      });
    }
};
