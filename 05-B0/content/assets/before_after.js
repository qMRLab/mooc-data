function beforeAfter() {
    $(".middle").on("mousedown", function () {
        $(document).on("mouseup",function() {
            $(document).unbind("mousemove")
             $(document).unbind("mouseup")
             oldX = undefined
             oldY = undefined
             })
        $(document).on("mousemove", function () {
            adjustSize()
        })
    })

    $(".before .modebar").css({'left':'0px', 'position':'absolute'})
}

window.fetch = new Proxy(window.fetch, {
    apply(fetch, that, args) {
        // Forward function call to the original fetch
        const result = fetch.apply(that, args);

        // Do whatever you want with the resulting Promise
        result.then((response) => {
            if (args[0] == '/_dash-update-component') {
                setTimeout(function() {beforeAfter()}, 500)
            }})
        return result
        }
    }
    )

var oldX;
var oldY;
function adjustSize() {
    if (event.screenX != 0) {
        var clientX = event.clientX
        var scrollLeft = $(document).scrollLeft()
        
        if (oldX && !((parseFloat($('.before').css("width"))+((clientX-scrollLeft)-oldX)) > $('.beforeAfter').width())) {
            var width = parseFloat($('.before').css("width"))+((clientX-scrollLeft)-oldX)
            $('.before').css("width", width)
            $('.middle').css("left", ($('.before').width()/$('.beforeAfter').width()*100)-2 + '%')
            
        }

        if (parseFloat($('.middle').css('left'))> $('.beforeAfter').width()) {
            $('.middle').css("left", '99.5%')
            $('.before').css("width", $('.beforeAfter').width())
        }
    }
    oldX = clientX-scrollLeft
    event.target.style.cursor = 'grab';
}

$(document).ready(function() {
    setTimeout(function() {beforeAfter()}, 1000)
})