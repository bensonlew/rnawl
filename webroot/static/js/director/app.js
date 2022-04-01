/*!
 * Author: Abdullah A Almsaeed
 * Date: 4 Jan 2014
 * Description:
 *      This file should be included in all pages
 !**/

/*
 * Global variables. If you change any of these vars, don't forget
 * to change the values in the less files!
 */
var left_side_width = 220; //Sidebar width in pixels

$(function() {
    "use strict";
    $('[data-toggle="popover"]').popover();
        $('.tree').treegrid();
        $('#jsonModal').on('show.bs.modal', function (event) {
            var code = $("#jsonModal pre>code.json")
            code.html("");
            var button = $(event.relatedTarget);
            $.post("json_data",{id:button.data("wid"),type:button.data("type")},function(data){
              code.html("<ul><li>" + data.replace(/\n/g,"\n</li><li>") +"\n</li></ul>");
              hljs.highlightBlock(code.get(0))
            });
        });
      $('button[data-toggle="dialog"]').click(function(e){
        BootstrapDialog.show({
        title:"WorkDir",
        message:$(this).data("content")
        })
      });
      $('button[data-toggle="dialog-cmd"]').click(function(e){
        BootstrapDialog.show({
        title:"Command",
        message:$(this).data("content")
        })
      });
      $('button[data-toggle="dialog-path"]').click(function(e){
        BootstrapDialog.show({
        title:"Path",
        message:$(this).data("content")
        })
      });
      $("pre>code").each(function(){
        $(this).html("<ul><li>" + $(this).html().replace(/\n/g,"\n</li><li>") +"\n</li></ul>");
      });
      hljs.initHighlightingOnLoad();
      if($("pre>code.python").length >0 && code_line>0){
          if($("pre>code.python>ul>li").length >= code_line){
              setTimeout(function(){
                  $("body,html").animate({
                  scrollTop:$("pre>code>ul>li:eq("+ code_line +")").offset().top - 100
                  }, 1000);
              },1000);
          }
      }

    $("#logout").click(function(e){
      window.location="logout";
    });
    $(".form_datetime").datetimepicker({format: 'yyyy-mm-dd hh:ii'});
    //Enable sidebar toggle
    $("[data-toggle='offcanvas']").click(function(e) {
        e.preventDefault();

        //If window is small enough, enable sidebar push menu
        if ($(window).width() <= 992) {
            $('.row-offcanvas').toggleClass('active');
            $('.left-side').removeClass("collapse-left");
            $(".right-side").removeClass("strech");
            $('.row-offcanvas').toggleClass("relative");
        } else {
            //Else, enable content streching
            $('.left-side').toggleClass("collapse-left");
            $(".right-side").toggleClass("strech");
        }
    });

    //Add hover support for touch devices
    $('.btn').bind('touchstart', function() {
        $(this).addClass('hover');
    }).bind('touchend', function() {
        $(this).removeClass('hover');
    });

  $('#main-menu').metisMenu();

  $(window).bind("load resize", function () {
    if ($(this).width() < 768) {
      $('div.sidebar-collapse').addClass('collapse')
    } else {
      $('div.sidebar-collapse').removeClass('collapse')
    }
  });


    function _fix() {
        //Get window height and the wrapper height
        var height = $(window).height() - $("body > .header").height() - ($("body > .footer").outerHeight() || 0);
        $(".wrapper").css("min-height", height + "px");
        var content = $(".wrapper").height();
        //If the wrapper height is greater than the window
        if (content > height)
            //then set sidebar height to the wrapper
            $(".left-side, html, body").css("min-height", content + "px");
        else {
            //Otherwise, set the sidebar to the height of the window
            $(".left-side, html, body").css("min-height", height + "px");
        }
    }
    //Fire upon load
    _fix();

});

/*
 * jQuery resize event - v1.1 - 3/14/2010
 * http://benalman.com/projects/jquery-resize-plugin/
 *
 * Copyright (c) 2010 "Cowboy" Ben Alman
 * Dual licensed under the MIT and GPL licenses.
 * http://benalman.com/about/license/
 */
(function($, h, c) {
    var a = $([]), e = $.resize = $.extend($.resize, {}), i, k = "setTimeout", j = "resize", d = j + "-special-event", b = "delay", f = "throttleWindow";
    e[b] = 250;
    e[f] = true;
    $.event.special[j] = {setup: function() {
            if (!e[f] && this[k]) {
                return false;
            }
            var l = $(this);
            a = a.add(l);
            $.data(this, d, {w: l.width(), h: l.height()});
            if (a.length === 1) {
                g();
            }
        }, teardown: function() {
            if (!e[f] && this[k]) {
                return false
            }
            var l = $(this);
            a = a.not(l);
            l.removeData(d);
            if (!a.length) {
                clearTimeout(i);
            }
        }, add: function(l) {
            if (!e[f] && this[k]) {
                return false
            }
            var n;
            function m(s, o, p) {
                var q = $(this), r = $.data(this, d);
                r.w = o !== c ? o : q.width();
                r.h = p !== c ? p : q.height();
                n.apply(this, arguments)
            }
            if ($.isFunction(l)) {
                n = l;
                return m
            } else {
                n = l.handler;
                l.handler = m
            }
        }};
    function g() {
        i = h[k](function() {
            a.each(function() {
                var n = $(this), m = n.width(), l = n.height(), o = $.data(this, d);
                if (m !== o.w || l !== o.h) {
                    n.trigger(j, [o.w = m, o.h = l])
                }
            });
            g()
        }, e[b])
    }}
)(jQuery, this);

/*
 * table search tool
 */
$(function() {
      data_type = $(".search-type-bar input[name='data_type']").val();
      set_bar($(".search-type-bar .dropdown-menu a[data-type='" + data_type +"']"), true);
      function set_bar(item, init){
          $(".search-type-bar button:first").html(item.text()+'<span class="caret"></span>');
          $(".search-type-bar input[name='data_type']").val(item.attr("data-type"));
          next = $(".search-type-bar .input-group-btn").next();
          data_type = item.attr("data-type");

          if(data_type =="run_id"||data_type=="path"||data_type=="task_id"||data_type=="api"){
              if(!next.is('input')){
                  next.replaceWith("<input type=\"text\" name=\"table_search\" class=\"form-control input-sm pull-right\" placeholder=\"Search\"/>");
                  return
              }else{
                $(".search-type-bar input[name='table_search']").val("")
              }
              if(init){
                    $(".search-type-bar input[name='table_search']").val(search_key)
              }
          }
          if(data_type =="add_time"||data_type =="run_time"||data_type =="end_time"){
              if(!next.is('div')){
                  next.replaceWith("<div class=\"row\"><div class=\"col-lg-6\"><input type=\"text\" readonly name=\"time_from\" class=\"form-control input-sm pull-right form_datetime col-lg-push-1\" placeholder=\"从...\"></div><div class=\"col-lg-6\"><input type=\"text\" readonly name=\"time_to\" class=\"form-control input-sm pull-right form_datetime col-lg-pull-1\" placeholder=\"到...\"></div></div>");
                  $(".form_datetime").datetimepicker({format: 'yyyy-mm-dd hh:ii'});

                  return
              }else{
                    $(".search-type-bar input[name='time_from']").val("")
                    $(".search-type-bar input[name='time_to']").val("")
              }
              if(init){
                    $(".search-type-bar input[name='time_from']").val(time_from)
                    $(".search-type-bar input[name='time_to']").val(time_to)
              }
          }

      }
      $(".search-type-bar .dropdown-menu a").click(function(e){
          e.preventDefault();
          set_bar($(this), false);
      })

  });


//  /*
//  *
//  * table control buttons
//  *
//  */
 $(function() {
    $("table.table .control-rerun").click(function(e){
        e.preventDefault();
        $("#rerun-form input[name='wid']").val($(this).data("wid"));
        $("#rerunModal .modal-footer .btn-primary").data("url", $(this).data("url"));
    });


  $("#rerunModal .modal-footer .btn-primary").click(function(e){
    e.preventDefault();
    btn = $(this);
    url = btn.data("url");
    btn.button('loading');
    form_data = $("#rerun-form").serialize();
    $.post(url, form_data,
       function(data){
         btn.button("reset");
         if(data.success){
            $("#rerunModal").modal('hide');
            d_type = BootstrapDialog.TYPE_SUCCESS;
            d_title = '操作成功';
            location.reload();
         }else{
            d_type = BootstrapDialog.TYPE_WARNING;
            d_title = '操作失败';
         }
         BootstrapDialog.show({
                    type: d_type,
                    title: d_title,
                    message: data.info,
                    buttons: [{
                        label: '关闭',
                        action: function(dialogRef){
                            dialogRef.close();
                        }
                    }]
         });
     }, "json");
  });
  function workflow_action(url, type, title, message, workflow_id, client){
    BootstrapDialog.show({
                    type: type,
                    title: title,
                    message: message,
                    buttons: [{
                        label: title,
                        action: function(dialogRef){
                            dialogRef.close()
                            $.post(url, {"client": client, "id": workflow_id},
                               function(data){
                                 if(data.success){
                                    d_type = BootstrapDialog.TYPE_SUCCESS;
                                    d_title = '操作成功';
                                 }else{
                                    d_type = BootstrapDialog.TYPE_WARNING;
                                    d_title = '操作失败';
                                 }
                                 BootstrapDialog.show({
                                            type: d_type,
                                            title: d_title,
                                            message: data.info,
                                            buttons: [{
                                                label: '关闭',
                                                action: function(dialogRef){
                                                    dialogRef.close();
                                                }
                                            }]
                                 });
                             }, "json");
                        }
                    },{
                        label: '关闭',
                        action: function(dialogRef){
                            dialogRef.close();
                        }
                    }]
         });
  }
  $("table.table .control-stop").click(function(e){
       type = BootstrapDialog.TYPE_DANGER
       msg = "确认停止流程的运行?此操作将不可撤销,操作成功后流程最多约需要15秒接收到指令，请刷新页面查看状态！"
       workflow_action("/pipeline/stop", type, "确认停止运行", msg,$(this).data("workflow"), $(this).data("client"))
  });
  $("table.table .control-continue").click(function(e){
       type = BootstrapDialog.TYPE_INFO
       msg = "确认继续暂停的流程?操作成功后流程最多约需要15秒接收到指令，请刷新页面查看状态！"
       workflow_action("/pipeline/stop_pause", type, "确认继续", msg,$(this).data("workflow"), $(this).data("client"))
  });
  $("table.table .control-pause").click(function(e){
       type = BootstrapDialog.TYPE_WARNING
       msg = "确认暂停流程的运行?须知暂停功能只是停止Tool任务的投递,暂停超过2小时,流程将自动退出运行！操作成功后流程最多约需要15秒接收到指令，请刷新页面查看状态！"
       workflow_action("/pipeline/pause", type, "确认暂停", msg,$(this).data("workflow"), $(this).data("client"))
  });


 });

/*
*
* log view page
*
*/
var current_code_index = null
var current_position = 0
var search_type = $('input:radio[name="logtype"]:checked').val();
function search_next_and_scroll(func,next){
    var last_node, first_node;
    if(current_code_index){
        current_code_index.removeClass("selected")
    }
    if(next){
        last_node = $("pre>code>ul>li:last")
        first_node = $("pre>code>ul>li:first")
        if(current_code_index==null){
            current_code_index = first_node
        }
        li_all = current_code_index.nextAll();
    }else{
        last_node = $("pre>code>ul>li:first")
        first_node = $("pre>code>ul>li:last")
        if(current_code_index==null){
            current_code_index = first_node
        }
        li_all = current_code_index.prevAll();
    }

    li_all.each(function(){
       if(func($(this).text())){
            $(this).addClass("selected");
            $("pre>code").animate({
                      scrollTop: 0
            }, 0);
            $("pre>code").animate({
                      scrollTop: $(this).position().top -200
            }, 0);

            current_code_index = $(this);
            return false
       }else{
        if($(this).is(last_node)){
            alert("已经搜索完毕！");
            current_code_index = first_node;
            current_position = 0;
        }
       }
    });
}

 $(function() {
    if(typeof(print_lines)!="undefined" && print_lines > 10){
        BootstrapDialog.show({
        type: BootstrapDialog.TYPE_WARNING,
        title:"警告",
        message:"该日志中有过多非日志格式输出，请注意调试完毕后关闭Print输出！<br/>必要的输出信息请使用self.logger方法输出，请不要输出过多不可阅读信息！"
        })
    };

        var regexp = new Array();
         regexp["error"]=/^\d{4}\-\d{2}\-\d{2} \d{2}:\d{2}:\d{2}\s+\S+\s+ERROR/i;
         regexp["warning"]=/^\d{4}\-\d{2}\-\d{2} \d{2}:\d{2}:\d{2}\s+\S+\s+WARNING/i;
         regexp["trace"]=/^Traceback/;
         regexp["print"]=/^Traceback/;
     var trace_start = false
     var key;
     function check(text){
        if(search_type!="print" && search_type!="key"){
            return text.match(regexp[search_type])
        }else if(search_type=="print"){
            if(text.match(/^\d{4}\-\d{2}\-\d{2} \d{2}:\d{2}/)){
                return false
            }else if(text.match(/^Traceback/)){
                trace_start = true
                return false
            }else{
                if(trace_start){
                    if(text.match(/^\w+:/)){
                        trace_start = false
                    }
                }else{
                    return true
                }
            }
        }else{
            if(text.indexOf(key)>=0){
                return true
            }else{
                return false
            }
        }
     }

    $("#findbar button.prev").click(function(){
        search_type = $('input:radio[name="logtype"]:checked').val();
        search_next_and_scroll(check, false);
    });
    $("#findbar button.next").click(function(){
        search_type = $('input:radio[name="logtype"]:checked').val();
        key = $("#findbar input[name='key_search']").val();
        if (search_type =="key" && key.length == 0){
            alert("请输入关键字！")
            return
        }
        search_next_and_scroll(check, true);
    });

 });

/*
*
* filter log bar
*
*/
 $(function() {
    $("#filterbar button.search").click(function(){
        filter_type = $(this).data("type");
        url = $(this).data("url")
        if(filter_type == "key"){
            key = $("#filterbar input[name='table_search']").val();
            if(key == ""){
                alert("请输入关键字！")
                return
            }else{
                $.post(url,{id:$(this).data("id"),key: key},function(data){
                  $("pre>code.log").html("<ul><li>" + data.replace(/\n/g,"\n</li><li>") +"\n</li></ul>");
                  hljs.highlightBlock($("pre>code.log").get(0));
                });
            }
        }else{
            time = $("#filterbar input[name='time_to']").val();
            area = $("#filterbar select[name='time_area']").val();
            $.post(url,{id:$(this).data("id"),time: time, area:area},function(data){
                  $("pre>code.log").html("<ul><li>" + data.replace(/\n/g,"\n</li><li>") +"\n</li></ul>");
                  hljs.highlightBlock($("pre>code.log").get(0));
            });
        }
    });
    $("#filterbar .dropdown-menu a").click(function(e){
          e.preventDefault();
          type = $(this).data("type");
          $("#filterbar button:first").text($(this).text());
          $("#filterbar button.search").data("type", type);
          next = $("#filterbar .input-group-btn").next();
          if(type=="key"){
            if(!next.is('input')){
                 next.replaceWith("<input type=\"text\" name=\"table_search\" class=\"form-control input-sm pull-right\" placeholder=\"Search\"/>");
                 return
            }
          }else{
            if(!next.is('div')){
                  next.replaceWith('<div class="row">'+
                          '<div class="col-lg-6">' +
                           '<input type="text" name="time_to" class="form-control input-sm pull-right form_datetime col-lg-push-1" placeholder="按时间筛选"></div>'+
                          '<div class="col-lg-6">' +
                          '<select name="time_area" class="form-control input-sm pull-right filter-time-area col-lg-pull-1">'+
                              '<option value=10 selected>前后10秒</option>'+
                              '<option value=20>前后20秒</option>'+
                              '<option value=30>前后30秒</option>'+
                              '<option value=60>前后1分钟</option>'+
                              '<option value=120>前后2分钟</option>'+
                              '<option value=180>前后3分钟</option>'+
                            '</select></div></div>');
                  $("#filterbar .form_datetime").datetimepicker({
                    format: 'yyyy-mm-dd hh:ii:ss',
                    minuteStep: 1,
                    todayHighlight:false,
                    initialDate:start_date,
                    startDate: start_date,
                    endDate: end_date
                  });
                  return
            }
          }

     })
 });

 /*
 *
 * api log view
 *
 */
 $(function(){
//    $("#apilog_search input[name='time_to']").datetimepicker({
//        format: 'yyyy-mm-dd hh:ii:ss',
//        minuteStep: 1,
//        todayHighlight:false,
//        initialDate:start_date,
//        startDate: start_date,
//        endDate: end_date
//    });
//    $("#apilog_search button.search").click(function(){
//        time = $("#apilog_search input[name='time_to']").val();
//        area = $("#apilog_search select[name='time_area']").val();
//        $.post("apilog",{workflow_id:$(this).data("wid"),time: time, area:area},function(data){
//              $("pre>code.apilog").html("<ul><li>" + data.replace(/\n/g,"\n</li><li>") +"\n</li></ul>");
//              hljs.highlightBlock($("pre>code.apilog").get(0));
//        });
//    });
       $('#ApiLogJsonModal').on('show.bs.modal', function (event) {
            var code = $("#ApiLogJsonModal pre>code.json")
            code.html("");
            var button = $(event.relatedTarget);
            data = JSON.stringify($.parseJSON(decodeURIComponent(button.data("json"))), null, 2)
            code.html("<ul><li>" + data.replace(/\n/g,"\n</li><li>") +"\n</li></ul>");
            hljs.highlightBlock(code.get(0))
       });
       $('#urlModal').on('show.bs.modal', function (event) {
            $('#urlModal div.modal-body ').find("a").remove();
            var button = $(event.relatedTarget);
            data = button.data("urlcode");
            $('#urlModal div.modal-body').html('<a style="word-wrap:break-word;" href="#" class="uri">'+ data +'</a>');
            $('#urlModal .modal-body').find("a").click(function(e){
                e.preventDefault();
                $("#ApiLogURLJsonModal").modal()
            });


            //alert();
       });
       $('#ApiLogURLJsonModal').on('show.bs.modal', function (event) {
            var code = $("#ApiLogURLJsonModal pre>code.json")
            code.html("");
            var text =$('#urlModal .modal-body').find("a").text();
            data = JSON.stringify($.parseJSON(decodeURIComponent(text.replace(/^sync_task_log=/,"")).replace(/\+/g, " ")), null, 2);
            code.html("<ul><li>" + data.replace(/\n/g,"\n</li><li>") +"\n</li></ul>");
            hljs.highlightBlock(code.get(0))
       });

 });