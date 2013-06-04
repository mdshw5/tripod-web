function wait_for_update(id) {
    $.ajax({ url: '/status/' + id,
             success: load_results, // if /status signals results ready, load them!
             complete: wait_for_update,  // if the wait_for_update poll times out, rerun
             timeout:  60000,
           });    
}

function load_results() {
    $("div#results").load("/results");
}
