/*! For license information please see 14.2439f61f.chunk.js.LICENSE.txt */
(this["webpackJsonpstreamlit-browser"]=this["webpackJsonpstreamlit-browser"]||[]).push([[14],{3663:function(t,e){},3664:function(t,e){},4321:function(t,e){},5755:function(t,e,i){"use strict";i.r(e),i.d(e,"default",(function(){return d}));var n=i(0),c=i(4301),r=i(4319),h=i(22),a=i(213),o=i(13),s=i.n(o)()("div",{target:"e1p558ko0"})((function(t){return{"& *":{fontFamily:t.theme.genericFonts.bodyFont,fontSize:"9.6px"},"& svg":{maxWidth:"100%"}}}),""),u=i(5);r.graphviz;var d=Object(a.a)((function(t){var e=t.width,i=t.element,r=t.height,a="graphviz-chart-".concat(i.elementId),o=0,d=0,g=function(){var t=d,n=o;return r?(t=e,n=r):i.useContainerWidth&&(t=e),{chartWidth:t,chartHeight:n}},f=function(){try{var t=Object(c.select)("#".concat(a)).graphviz().zoom(!1).fit(!0).scale(1).renderDot(i.spec).on("end",(function(){var t=Object(c.select)("#".concat(a," > svg")).node();t&&(o=t.getBBox().height,d=t.getBBox().width)})),e=g(),n=e.chartHeight,r=e.chartWidth;n>0&&t.height(n),r>0&&t.width(r)}catch(s){Object(h.b)(s)}};Object(n.useEffect)((function(){f()}));var v=g(),p=v.chartWidth?v.chartWidth:e,b=v.chartHeight?v.chartHeight:r;return Object(u.jsx)(s,{className:"graphviz stGraphVizChart",id:a,style:{width:p,height:b}})}))}}]);