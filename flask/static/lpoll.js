function wait_for_update(pid) {
    $.ajax({ url: '/status/' + pid,
             success: $("div#results").load("/results");, // if /status signals results ready, load them!
             complete: wait_for_update,  // if the wait_for_update poll times out, rerun
             timeout:  60000,
           });    
}

