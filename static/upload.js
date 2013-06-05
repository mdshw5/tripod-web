$(document).ready(function(){
    $("#submitted").hide();
    $(".extras").hide();
    $("#options").submit(function(){
	$("#submitted").toggle();
    });
    $("#show_hide").click(function(){
	$(".extras").slideToggle();
    });
});

