function wait_for_update(id) {
    $.ajax({ url: '/status/' + id,
	     success: redirect,
	     dataType: "text",
             complete: wait_for_update(id),
             timeout:  60000 });    
}

function redirect() {
    window.location.replace('/results');
}