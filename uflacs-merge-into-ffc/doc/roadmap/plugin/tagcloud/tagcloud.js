// tagcloud
(function(){
[].forEach.call( document.querySelectorAll('[tagcloud]'), function(cloud){
    cloud.innerHTML = '<span>' + cloud.innerHTML.split(/\n/).join('</span> <span>') + '</span>';
    [].forEach.call( cloud.querySelectorAll('span'), function(elem){
        var prctnge = Math.random() * 150 + 50;
        if (cloud.hasAttribute('large')) {
            prctnge = prctnge * 1.2;
        }
        elem.style.fontSize = prctnge + '%';
        if (cloud.hasAttribute('bw')) {
            var col = Math.round(Math.random() * 155 + 100);
            elem.style.color = 'rgb('+ col  +',' + col + ',' + col + ')'
        } else {
            elem.style.color = 'hsl('+ Math.random()*360 +', 40%, 50%)'
        }
        elem.classList.add('clouditem')
    });
});
}
)();