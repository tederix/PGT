(this["webpackJsonpstreamlit-browser"]=this["webpackJsonpstreamlit-browser"]||[]).push([[37],{5769:function(e,t,o){"use strict";o.r(t),o.d(t,"default",(function(){return c}));var r=o(6),n=o(9),a=o(7),i=o(8),l=o(0),u=o.n(l),s=o(211),m=o(172),p=o(5),c=function(e){Object(a.a)(o,e);var t=Object(i.a)(o);function o(){var e;Object(r.a)(this,o);for(var n=arguments.length,a=new Array(n),i=0;i<n;i++)a[i]=arguments[i];return(e=t.call.apply(t,[this].concat(a))).formClearHelper=new s.b,e.state={value:e.initialValue},e.commitWidgetValue=function(t){e.props.widgetMgr.setStringValue(e.props.element,e.state.value,t)},e.onFormCleared=function(){e.setState({value:e.props.element.default},(function(){return e.commitWidgetValue({fromUi:!0})}))},e.onColorClose=function(t){e.setState({value:t},(function(){return e.commitWidgetValue({fromUi:!0})}))},e.render=function(){var t=e.props,o=t.element,r=t.width,n=t.disabled,a=t.widgetMgr,i=e.state.value;return e.formClearHelper.manageFormClearListener(a,o.formId,e.onFormCleared),Object(p.jsx)(m.a,{label:o.label,help:o.help,onChange:e.onColorClose,disabled:n,width:r,value:i})},e}return Object(n.a)(o,[{key:"initialValue",get:function(){var e=this.props.widgetMgr.getStringValue(this.props.element);return void 0!==e?e:this.props.element.default}},{key:"componentDidMount",value:function(){this.props.element.setValue?this.updateFromProtobuf():this.commitWidgetValue({fromUi:!1})}},{key:"componentDidUpdate",value:function(){this.maybeUpdateFromProtobuf()}},{key:"componentWillUnmount",value:function(){this.formClearHelper.disconnect()}},{key:"maybeUpdateFromProtobuf",value:function(){this.props.element.setValue&&this.updateFromProtobuf()}},{key:"updateFromProtobuf",value:function(){var e=this,t=this.props.element.value;this.props.element.setValue=!1,this.setState({value:t},(function(){e.commitWidgetValue({fromUi:!1})}))}}]),o}(u.a.PureComponent)}}]);