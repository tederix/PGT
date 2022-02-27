/*! For license information please see 10.28cc38f0.chunk.js.LICENSE.txt */
(this["webpackJsonpstreamlit-browser"]=this["webpackJsonpstreamlit-browser"]||[]).push([[10],{3485:function(t,e){},3490:function(t,e){},3645:function(t,e){},3648:function(t,e){},3801:function(t,e){},3817:function(t,e){},3819:function(t,e){},3846:function(t,e){},4247:function(t,e){},5748:function(t,e,n){"use strict";n.r(e),n.d(e,"default",(function(){return K}));var i=n(23),a=n(6),r=n(9),o=n(7),c=n(8),s=n(10),h=n(0),u=n(4291),l=n.n(u),p=n(3743),b=n.n(p),j=n(3744),d=n(3197),f=n(4282),w=n(3583),m=n(3578),O=n(3341),v=n(4220),g=n(4233),x=n(3198),k=n(213),S=n(11),y=n.n(S),T=n(17),V=n(137),L=n(56),C=n(5734),M=n(216),E=n.n(M),F=n(38),N=function(t){Object(o.a)(n,t);var e=Object(c.a)(n);function n(){return Object(a.a)(this,n),e.apply(this,arguments)}return n}(Object(C.a)(Error)),D=function(t){Object(o.a)(n,t);var e=Object(c.a)(n);function n(){return Object(a.a)(this,n),e.apply(this,arguments)}return n}(Object(C.a)(Error)),J=function(){function t(){Object(a.a)(this,t)}return Object(r.a)(t,null,[{key:"get",value:function(){var e=Object(T.a)(y.a.mark((function e(){var n,i,a;return y.a.wrap((function(e){for(;;)switch(e.prev=e.next){case 0:if(n=F.a.current,i=n.commandLine,a=n.userMapboxToken,t.token&&t.commandLine===i.toLowerCase()){e.next=10;break}if(""===a){e.next=6;break}t.token=a,e.next=9;break;case 6:return e.next=8,this.fetchToken("https://data.streamlit.io/tokens.json","mapbox");case 8:t.token=e.sent;case 9:t.commandLine=i.toLowerCase();case 10:return e.abrupt("return",t.token);case 11:case"end":return e.stop()}}),e,this)})));return function(){return e.apply(this,arguments)}}()},{key:"fetchToken",value:function(){var t=Object(T.a)(y.a.mark((function t(e,n){var i,a;return y.a.wrap((function(t){for(;;)switch(t.prev=t.next){case 0:return t.prev=0,t.next=3,E.a.get(e);case 3:if(i=t.sent,null!=(a=i.data[n])&&""!==a){t.next=7;break}throw new Error('Missing token "'.concat(n,'"'));case 7:return t.abrupt("return",a);case 10:throw t.prev=10,t.t0=t.catch(0),new D("".concat(t.t0.message," (").concat(e,")"));case 13:case"end":return t.stop()}}),t,null,[[0,10]])})));return function(e,n){return t.apply(this,arguments)}}()}]),t}();J.token=void 0,J.commandLine=void 0,J.isRunningLocal=function(){var t=window.location.hostname;return"localhost"===t||"127.0.0.1"===t};var P=n(126),z=n.n(P),A=n(138),I=n(5),q=function(t){var e=t.error,n=t.width,i=t.deltaType;return e instanceof N?Object(I.jsx)(A.a,{width:n,name:"No Mapbox token provided",message:Object(I.jsxs)(I.Fragment,{children:[Object(I.jsxs)("p",{children:["To use ",Object(I.jsxs)("code",{children:["st.",i]})," or ",Object(I.jsx)("code",{children:"st.map"})," you need to set up a Mapbox access token."]}),Object(I.jsxs)("p",{children:["To get a token, create an account at"," ",Object(I.jsx)("a",{href:"https://mapbox.com",children:"https://mapbox.com"}),". It's free for moderate usage levels!"]}),Object(I.jsxs)("p",{children:["Once you have a token, just set it using the Streamlit config option ",Object(I.jsx)("code",{children:"mapbox.token"})," and don't forget to restart your Streamlit server at this point if it's still running, then reload this tab."]}),Object(I.jsxs)("p",{children:["See"," ",Object(I.jsx)("a",{href:"https://docs.streamlit.io/library/advanced-features/configuration#view-all-configuration-options",children:"our documentation"})," ","for more info on how to set config options."]})]})}):e instanceof D?Object(I.jsx)(A.a,{width:n,name:"Error fetching Streamlit Mapbox token",message:Object(I.jsxs)(I.Fragment,{children:[Object(I.jsx)("p",{children:"This app requires an internet connection."}),Object(I.jsx)("p",{children:"Please check your connection and try again."}),Object(I.jsxs)("p",{children:["If you think this is a bug, please file bug report"," ",Object(I.jsx)("a",{href:"https://github.com/streamlit/streamlit/issues/new/choose",children:"here"}),"."]})]})}):Object(I.jsx)(A.a,{width:n,name:"Error fetching Streamlit Mapbox token",message:e.message})},G=function(t){return function(e){var n=function(n){Object(o.a)(r,n);var i=Object(c.a)(r);function r(n){var o;return Object(a.a)(this,r),(o=i.call(this,n)).initMapboxToken=Object(T.a)(y.a.mark((function t(){var e;return y.a.wrap((function(t){for(;;)switch(t.prev=t.next){case 0:return t.prev=0,t.next=3,J.get();case 3:e=t.sent,o.setState({mapboxToken:e,isFetching:!1}),t.next=10;break;case 7:t.prev=7,t.t0=t.catch(0),o.setState({mapboxTokenError:t.t0,isFetching:!1});case 10:case"end":return t.stop()}}),t,null,[[0,7]])}))),o.render=function(){var n=o.state,i=n.mapboxToken,a=n.mapboxTokenError,r=n.isFetching,c=o.props.width;return a?Object(I.jsx)(q,{width:c,error:a,deltaType:t}):r?Object(I.jsx)(V.a,{body:"Loading...",kind:L.a.INFO,width:c}):Object(I.jsx)(e,Object(s.a)({mapboxToken:i},o.props))},o.state={isFetching:!0,mapboxToken:void 0,mapboxTokenError:void 0},o.initMapboxToken(),o}return r}(h.PureComponent);return n.displayName="withMapboxToken(".concat(e.displayName||e.name,")"),z()(n,e)}},R=n(13),W=n.n(R)()("div",{target:"e1q4dr930"})((function(t){var e=t.width;return{position:"relative",height:t.height,width:e}}),""),_=(n(3792),{classes:Object(s.a)(Object(s.a)(Object(s.a)(Object(s.a)({},d),m),w),O)});Object(x.registerLoaders)([v.CSVLoader,g.GLTFLoader]);var B=new f.JSONConverter({configuration:_}),H=function(t){Object(o.a)(n,t);var e=Object(c.a)(n);function n(){var t;Object(a.a)(this,n);for(var i=arguments.length,r=new Array(i),o=0;o<i;o++)r[o]=arguments[o];return(t=e.call.apply(e,[this].concat(r))).state={viewState:{bearing:0,pitch:0,zoom:11},initialized:!1,initialViewState:{}},t.componentDidMount=function(){t.setState({initialized:!0})},t.createTooltip=function(e){var n=t.props.element;if(!e||!e.object||!n.tooltip)return!1;var i=JSON.parse(n.tooltip);return i.html?i.html=t.interpolate(e,i.html):i.text=t.interpolate(e,i.text),i},t.interpolate=function(t,e){var n=e.match(/{(.*?)}/g);return n&&n.forEach((function(n){var i=n.substring(1,n.length-1);t.object.hasOwnProperty(i)&&(e=e.replace(n,t.object[i]))})),e},t.onViewStateChange=function(e){var n=e.viewState;t.setState({viewState:n})},t}return Object(r.a)(n,[{key:"render",value:function(){var t=n.getDeckObject(this.props),e=this.state.viewState;return Object(I.jsx)(W,{className:"stDeckGlJsonChart",width:t.initialViewState.width,height:t.initialViewState.height,children:Object(I.jsx)(l.a,{viewState:e,onViewStateChange:this.onViewStateChange,height:t.initialViewState.height,width:t.initialViewState.width,layers:this.state.initialized?t.layers:[],getTooltip:this.createTooltip,controller:!0,children:Object(I.jsx)(j.StaticMap,{height:t.initialViewState.height,width:t.initialViewState.width,mapStyle:t.mapStyle&&("string"===typeof t.mapStyle?t.mapStyle:t.mapStyle[0]),mapboxApiAccessToken:this.props.mapboxToken})})})}}],[{key:"getDerivedStateFromProps",value:function(t,e){var a=n.getDeckObject(t);if(!b()(a.initialViewState,e.initialViewState)){var r=Object.keys(a.initialViewState).reduce((function(t,n){return a.initialViewState[n]===e.initialViewState[n]?t:Object(s.a)(Object(s.a)({},t),{},Object(i.a)({},n,a.initialViewState[n]))}),{});return{viewState:Object(s.a)(Object(s.a)({},e.viewState),r),initialViewState:a.initialViewState}}return null}}]),n}(h.PureComponent);H.getDeckObject=function(t){var e=t.element,n=t.width,i=t.height,a=JSON.parse(e.json);return i?(a.initialViewState.height=i,a.initialViewState.width=n):(a.initialViewState.height||(a.initialViewState.height=500),e.useContainerWidth&&(a.initialViewState.width=n)),delete a.views,B.convert(a)};var K=G("st.pydeck_chart")(Object(k.a)(H))}}]);