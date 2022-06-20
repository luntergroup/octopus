"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[838],{3905:function(e,r,n){n.d(r,{Zo:function(){return p},kt:function(){return m}});var t=n(7294);function o(e,r,n){return r in e?Object.defineProperty(e,r,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[r]=n,e}function a(e,r){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var t=Object.getOwnPropertySymbols(e);r&&(t=t.filter((function(r){return Object.getOwnPropertyDescriptor(e,r).enumerable}))),n.push.apply(n,t)}return n}function i(e){for(var r=1;r<arguments.length;r++){var n=null!=arguments[r]?arguments[r]:{};r%2?a(Object(n),!0).forEach((function(r){o(e,r,n[r])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):a(Object(n)).forEach((function(r){Object.defineProperty(e,r,Object.getOwnPropertyDescriptor(n,r))}))}return e}function s(e,r){if(null==e)return{};var n,t,o=function(e,r){if(null==e)return{};var n,t,o={},a=Object.keys(e);for(t=0;t<a.length;t++)n=a[t],r.indexOf(n)>=0||(o[n]=e[n]);return o}(e,r);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);for(t=0;t<a.length;t++)n=a[t],r.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(o[n]=e[n])}return o}var c=t.createContext({}),l=function(e){var r=t.useContext(c),n=r;return e&&(n="function"==typeof e?e(r):i(i({},r),e)),n},p=function(e){var r=l(e.components);return t.createElement(c.Provider,{value:r},e.children)},u={inlineCode:"code",wrapper:function(e){var r=e.children;return t.createElement(t.Fragment,{},r)}},d=t.forwardRef((function(e,r){var n=e.components,o=e.mdxType,a=e.originalType,c=e.parentName,p=s(e,["components","mdxType","originalType","parentName"]),d=l(n),m=o,f=d["".concat(c,".").concat(m)]||d[m]||u[m]||a;return n?t.createElement(f,i(i({ref:r},p),{},{components:n})):t.createElement(f,i({ref:r},p))}));function m(e,r){var n=arguments,o=r&&r.mdxType;if("string"==typeof e||o){var a=n.length,i=new Array(a);i[0]=d;var s={};for(var c in r)hasOwnProperty.call(r,c)&&(s[c]=r[c]);s.originalType=e,s.mdxType="string"==typeof e?e:o,i[1]=s;for(var l=2;l<a;l++)i[l]=n[l];return t.createElement.apply(null,i)}return t.createElement.apply(null,n)}d.displayName="MDXCreateElement"},8278:function(e,r,n){n.r(r),n.d(r,{default:function(){return u},frontMatter:function(){return s},metadata:function(){return c},toc:function(){return l}});var t=n(7462),o=n(3366),a=(n(7294),n(3905)),i=["components"],s={id:"errorModels",title:"Error Models"},c={unversionedId:"guides/errorModels",id:"guides/errorModels",isDocsHomePage:!1,title:"Error Models",description:"Octopus accounts for SNV and indel sequencing errors with a context aware error model. The parameterisation of this model is conditional on the library preparations and sequencing technology used, and can have consequences on calling accuracy, particular for indel errors in tandem repeat regions. Octopus comes packaged with parameter sets for several common library preparation and sequencing combinations, and also allows custom sequence error models to be used.",source:"@site/docs/guides/errorModels.md",sourceDirName:"guides",slug:"/guides/errorModels",permalink:"/octopus/docs/guides/errorModels",editUrl:"https://github.com/${organizationName}/${projectName}/edit/${branch}/website/docs/guides/errorModels.md",version:"current",frontMatter:{id:"errorModels",title:"Error Models"},sidebar:"docs",previous:{title:"Cell",permalink:"/octopus/docs/guides/models/cell"},next:{title:"Phasing",permalink:"/octopus/docs/guides/phasing"}},l=[],p={toc:l};function u(e){var r=e.components,n=(0,o.Z)(e,i);return(0,a.kt)("wrapper",(0,t.Z)({},p,n,{components:r,mdxType:"MDXLayout"}),(0,a.kt)("p",null,"Octopus accounts for SNV and indel sequencing errors with a context aware error model. The parameterisation of this model is conditional on the library preparations and sequencing technology used, and can have consequences on calling accuracy, particular for indel errors in tandem repeat regions. Octopus comes packaged with parameter sets for several common library preparation and sequencing combinations, and also allows custom sequence error models to be used."),(0,a.kt)("p",null,"Built-in error models are selected using the ",(0,a.kt)("inlineCode",{parentName:"p"},"--sequence-error-model")," option, which accepts inputs of the form ",(0,a.kt)("inlineCode",{parentName:"p"},"[library preparation]<.sequencer>"),". library preparation is selected from: ",(0,a.kt)("inlineCode",{parentName:"p"},"PCR"),", ",(0,a.kt)("inlineCode",{parentName:"p"},"PCR-FREE"),", or ",(0,a.kt)("inlineCode",{parentName:"p"},"10X"),". sequencer is selected from: ",(0,a.kt)("inlineCode",{parentName:"p"},"HISEQ-2000"),", ",(0,a.kt)("inlineCode",{parentName:"p"},"HISEQ-2500"),", ",(0,a.kt)("inlineCode",{parentName:"p"},"HISEQ-4000"),", ",(0,a.kt)("inlineCode",{parentName:"p"},"X10"),", ",(0,a.kt)("inlineCode",{parentName:"p"},"NOVASEQ"),", ",(0,a.kt)("inlineCode",{parentName:"p"},"BGISEQ-5000"),". For example, ",(0,a.kt)("inlineCode",{parentName:"p"},"PCR.NOVASEQ")," would select the sequence error model parametrised for a ",(0,a.kt)("inlineCode",{parentName:"p"},"PCR")," library preparation and a ",(0,a.kt)("inlineCode",{parentName:"p"},"NOVASEQ")," sequencer. If no sequencer is provided then the default is used (see ",(0,a.kt)("inlineCode",{parentName:"p"},"octopus --help"),")."),(0,a.kt)("p",null,"Custom error models can be used by providing a path to a valid Octopus error model file. These can be produced using the ",(0,a.kt)("a",{parentName:"p",href:"https://github.com/luntergroup/octopus/blob/develop/scripts/profiler.py"},(0,a.kt)("inlineCode",{parentName:"a"},"profiler.py"))," Python script in the scripts top level directory. The script creates error model files given the output of the ",(0,a.kt)("inlineCode",{parentName:"p"},"--data-profile")," command line option."))}u.isMDXComponent=!0}}]);