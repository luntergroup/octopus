"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[217],{3905:function(e,t,n){n.d(t,{Zo:function(){return c},kt:function(){return m}});var a=n(7294);function i(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function r(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);t&&(a=a.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,a)}return n}function o(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?r(Object(n),!0).forEach((function(t){i(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):r(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function l(e,t){if(null==e)return{};var n,a,i=function(e,t){if(null==e)return{};var n,a,i={},r=Object.keys(e);for(a=0;a<r.length;a++)n=r[a],t.indexOf(n)>=0||(i[n]=e[n]);return i}(e,t);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);for(a=0;a<r.length;a++)n=r[a],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(i[n]=e[n])}return i}var s=a.createContext({}),p=function(e){var t=a.useContext(s),n=t;return e&&(n="function"==typeof e?e(t):o(o({},t),e)),n},c=function(e){var t=p(e.components);return a.createElement(s.Provider,{value:t},e.children)},u={inlineCode:"code",wrapper:function(e){var t=e.children;return a.createElement(a.Fragment,{},t)}},d=a.forwardRef((function(e,t){var n=e.components,i=e.mdxType,r=e.originalType,s=e.parentName,c=l(e,["components","mdxType","originalType","parentName"]),d=p(n),m=i,h=d["".concat(s,".").concat(m)]||d[m]||u[m]||r;return n?a.createElement(h,o(o({ref:t},c),{},{components:n})):a.createElement(h,o({ref:t},c))}));function m(e,t){var n=arguments,i=t&&t.mdxType;if("string"==typeof e||i){var r=n.length,o=new Array(r);o[0]=d;var l={};for(var s in t)hasOwnProperty.call(t,s)&&(l[s]=t[s]);l.originalType=e,l.mdxType="string"==typeof e?e:i,o[1]=l;for(var p=2;p<r;p++)o[p]=n[p];return a.createElement.apply(null,o)}return a.createElement.apply(null,n)}d.displayName="MDXCreateElement"},872:function(e,t,n){n.r(t),n.d(t,{default:function(){return u},frontMatter:function(){return l},metadata:function(){return s},toc:function(){return p}});var a=n(7462),i=n(3366),r=(n(7294),n(3905)),o=["components"],l={id:"installation",title:"Installation"},s={unversionedId:"installation",id:"installation",isDocsHomePage:!1,title:"Installation",description:"Octopus can be built and installed on most Unix based systems (e.g. Linux and MacOS). Windows has not been tested.",source:"@site/docs/installation.md",sourceDirName:".",slug:"/installation",permalink:"/octopus/docs/installation",editUrl:"https://github.com/${organizationName}/${projectName}/edit/${branch}/website/docs/installation.md",version:"current",frontMatter:{id:"installation",title:"Installation"},sidebar:"docs",previous:{title:"Introduction",permalink:"/octopus/docs/introduction"},next:{title:"Overview",permalink:"/octopus/docs/guides/overview"}},p=[{value:"Requirements",id:"requirements",children:[]},{value:"Python",id:"python",children:[{value:"Setting the build architecture",id:"setting-the-build-architecture",children:[]}]},{value:"CMake",id:"cmake",children:[]},{value:"Docker",id:"docker",children:[]},{value:"Singularity",id:"singularity",children:[]},{value:"Conda",id:"conda",children:[]}],c={toc:p};function u(e){var t=e.components,n=(0,i.Z)(e,o);return(0,r.kt)("wrapper",(0,a.Z)({},c,n,{components:t,mdxType:"MDXLayout"}),(0,r.kt)("p",null,"Octopus can be built and installed on most Unix based systems (e.g. Linux and MacOS). Windows has not been tested."),(0,r.kt)("p",null,"The recommend way to install Octopus for most users is:"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-shell"},"$ git clone -b master https://github.com/luntergroup/octopus.git\n$ octopus/scripts/install.py --dependencies --forests\n")),(0,r.kt)("p",null,"You can then optionally add ",(0,r.kt)("inlineCode",{parentName:"p"},"octopus")," to your ",(0,r.kt)("inlineCode",{parentName:"p"},"PATH"),":"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-shell"},"$ echo 'export PATH='$(pwd)'/octopus/bin:$PATH' >> ~/.bash_profile\n$ source ~/.bash_profile\n")),(0,r.kt)("p",null,"Then check the installation was successful:"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-shell"},"$ octopus --version\n")),(0,r.kt)("h2",{id:"requirements"},"Requirements"),(0,r.kt)("ul",null,(0,r.kt)("li",{parentName:"ul"},"A ",(0,r.kt)("a",{parentName:"li",href:"https://isocpp.org/wiki/faq/cpp14"},"C++14")," compiler and compatibility standard library. Either ",(0,r.kt)("a",{parentName:"li",href:"https://gcc.gnu.org"},"GCC")," (version >= 9.3) or ",(0,r.kt)("a",{parentName:"li",href:"https://clang.llvm.org"},"Clang")," (version >= 11.0) are recommended."),(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("a",{parentName:"li",href:"https://git-scm.com"},"Git")," version >= 2.5"),(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("a",{parentName:"li",href:"https://www.boost.org"},"Boost")," version >= 1.65"),(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("a",{parentName:"li",href:"https://github.com/samtools/htslib"},"htslib")," version >= 1.4; version != 1.12"),(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("a",{parentName:"li",href:"https://gmplib.org"},"GMP")," version >= 5.1.0"),(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("a",{parentName:"li",href:"https://cmake.org"},"CMake")," version >= 3.9"),(0,r.kt)("li",{parentName:"ul"},"Optional:",(0,r.kt)("ul",{parentName:"li"},(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("a",{parentName:"li",href:"https://www.python.org"},"Python")," version >= 3 plus the ",(0,r.kt)("a",{parentName:"li",href:"https://pypi.org/project/distro/"},"distro")," package")))),(0,r.kt)("div",{className:"admonition admonition-important alert alert--info"},(0,r.kt)("div",{parentName:"div",className:"admonition-heading"},(0,r.kt)("h5",{parentName:"div"},(0,r.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,r.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"14",height:"16",viewBox:"0 0 14 16"},(0,r.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M7 2.3c3.14 0 5.7 2.56 5.7 5.7s-2.56 5.7-5.7 5.7A5.71 5.71 0 0 1 1.3 8c0-3.14 2.56-5.7 5.7-5.7zM7 1C3.14 1 0 4.14 0 8s3.14 7 7 7 7-3.14 7-7-3.14-7-7-7zm1 3H6v5h2V4zm0 6H6v2h2v-2z"}))),"important")),(0,r.kt)("div",{parentName:"div",className:"admonition-content"},(0,r.kt)("p",{parentName:"div"},"Octopus uses ",(0,r.kt)("a",{parentName:"p",href:"https://en.wikipedia.org/wiki/SIMD"},"SIMD")," instructions for performance reasons. The instruction set used (minimum ",(0,r.kt)("a",{parentName:"p",href:"https://en.wikipedia.org/wiki/SSE2"},"SSE2"),") is built statically, so if you compile with ",(0,r.kt)("a",{parentName:"p",href:"https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#Advanced_Vector_Extensions_2"},"AVX2"),", you won't be able to use the resulting binary on machines that doesn't support AVX2."))),(0,r.kt)("h2",{id:"python"},"Python"),(0,r.kt)("p",null,"First clone the git repository in your preferred directory:"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-shell"},"$ git clone -b master https://github.com/luntergroup/octopus.git && cd octopus\n")),(0,r.kt)("p",null,"The easiest way to install Octopus from source is with the Python3 installer script. To see the options available to this script run ",(0,r.kt)("inlineCode",{parentName:"p"},"scripts/install.py --help"),"."),(0,r.kt)("p",null,"If all the requirements are accessible on your ",(0,r.kt)("inlineCode",{parentName:"p"},"PATH")," then simply run"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-shell"},"$ scripts/install.py\n")),(0,r.kt)("p",null,"otherwise you can specify paths to each dependency, for example, to set the compiler you'd use "),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-shell"},"$ scripts/install.py --cxx_compiler /path/to/cpp/compiler\n")),(0,r.kt)("p",null,"By default, this installs to ",(0,r.kt)("inlineCode",{parentName:"p"},"/bin")," relative to where octopus is installed. To install to a different location (e.g. ",(0,r.kt)("inlineCode",{parentName:"p"},"/usr/local/bin"),") use:"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-shell"},"$ scripts/install.py --prefix /user/local/bin\n")),(0,r.kt)("p",null,"You can also request all dependencies to be installed locally:"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-shell"},"$ scripts/install.py --dependencies\n")),(0,r.kt)("div",{className:"admonition admonition-tip alert alert--success"},(0,r.kt)("div",{parentName:"div",className:"admonition-heading"},(0,r.kt)("h5",{parentName:"div"},(0,r.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,r.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"12",height:"16",viewBox:"0 0 12 16"},(0,r.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M6.5 0C3.48 0 1 2.19 1 5c0 .92.55 2.25 1 3 1.34 2.25 1.78 2.78 2 4v1h5v-1c.22-1.22.66-1.75 2-4 .45-.75 1-2.08 1-3 0-2.81-2.48-5-5.5-5zm3.64 7.48c-.25.44-.47.8-.67 1.11-.86 1.41-1.25 2.06-1.45 3.23-.02.05-.02.11-.02.17H5c0-.06 0-.13-.02-.17-.2-1.17-.59-1.83-1.45-3.23-.2-.31-.42-.67-.67-1.11C2.44 6.78 2 5.65 2 5c0-2.2 2.02-4 4.5-4 1.22 0 2.36.42 3.22 1.19C10.55 2.94 11 3.94 11 5c0 .66-.44 1.78-.86 2.48zM4 14h5c-.23 1.14-1.3 2-2.5 2s-2.27-.86-2.5-2z"}))),"tip")),(0,r.kt)("div",{parentName:"div",className:"admonition-content"},(0,r.kt)("p",{parentName:"div"},"If a build isn't working after an update then try adding ",(0,r.kt)("inlineCode",{parentName:"p"},"--clean")," to the install command."))),(0,r.kt)("h3",{id:"setting-the-build-architecture"},"Setting the build architecture"),(0,r.kt)("p",null,"By default, the binary is optimised for the build machine architecture. If you need to run Octopus on another machine with a different architecture then use the ",(0,r.kt)("inlineCode",{parentName:"p"},"--architecture")," option:"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-shell"},"$ scripts/install.py --architecture haswell\n")),(0,r.kt)("p",null,"This is passed to the ",(0,r.kt)("a",{parentName:"p",href:"https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html"},"-march")," compiler option. "),(0,r.kt)("h2",{id:"cmake"},"CMake"),(0,r.kt)("p",null,"If Python3 isn't available, Octopus can be installed directly with ",(0,r.kt)("a",{parentName:"p",href:"https://cmake.org"},"CMake"),":"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-shell"},"$ git clone -b master https://github.com/luntergroup/octopus.git\n$ cd octopus/build\n$ cmake .. && make install\n")),(0,r.kt)("p",null,"CMake will try to find a suitable compiler on your system, if you'd like you use a specific compiler use the ",(0,r.kt)("inlineCode",{parentName:"p"},"-D")," option, for example:"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-shell"},"$ cmake -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ ..\n")),(0,r.kt)("h2",{id:"docker"},"Docker"),(0,r.kt)("p",null,"Pre-built Docker images are available on ",(0,r.kt)("a",{parentName:"p",href:"https://hub.docker.com/r/dancooke/octopus"},"DockerHub"),":"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-shell"},"$ docker pull dancooke/octopus\n$ docker run dancooke/octopus -h\n")),(0,r.kt)("div",{className:"admonition admonition-important alert alert--info"},(0,r.kt)("div",{parentName:"div",className:"admonition-heading"},(0,r.kt)("h5",{parentName:"div"},(0,r.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,r.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"14",height:"16",viewBox:"0 0 14 16"},(0,r.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M7 2.3c3.14 0 5.7 2.56 5.7 5.7s-2.56 5.7-5.7 5.7A5.71 5.71 0 0 1 1.3 8c0-3.14 2.56-5.7 5.7-5.7zM7 1C3.14 1 0 4.14 0 8s3.14 7 7 7 7-3.14 7-7-3.14-7-7-7zm1 3H6v5h2V4zm0 6H6v2h2v-2z"}))),"important")),(0,r.kt)("div",{parentName:"div",className:"admonition-content"},(0,r.kt)("p",{parentName:"div"},"The Octopus images on DockerHub are built to a ",(0,r.kt)("a",{parentName:"p",href:"https://en.wikipedia.org/wiki/Haswell_(microarchitecture)"},"Haswell")," architecture. This means that they will only work on  Haswell (with AVX2) or newer machines."))),(0,r.kt)("p",null,"You can also build a new image from the Dockerfile:"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-shell"},"$ git clone https://github.com/luntergroup/octopus.git && cd octopus\n$ docker build -t octopus .\n$ docker run octopus -h\n")),(0,r.kt)("p",null,"This is especially useful if you need to build to a specific architecture:"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-shell"},"$ docker build -t octopus --build-args architecture=sandybridge .\n")),(0,r.kt)("h2",{id:"singularity"},"Singularity"),(0,r.kt)("p",null,"To build a ",(0,r.kt)("a",{parentName:"p",href:"https://singularity.hpcng.org"},"Singularity")," container directly from the DockerHub images use"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-shell"},"$ singularity build octopus.sif docker://dancooke/octopus\n")),(0,r.kt)("h2",{id:"conda"},"Conda"),(0,r.kt)("p",null,"Octopus is available ",(0,r.kt)("a",{parentName:"p",href:"https://anaconda.org/bioconda/octopus"},"pre-built for Linux")," as part of ",(0,r.kt)("a",{parentName:"p",href:"https://bioconda.github.io/"},"Bioconda"),":"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-shell"},"$ conda install -c bioconda octopus\n")),(0,r.kt)("div",{className:"admonition admonition-important alert alert--info"},(0,r.kt)("div",{parentName:"div",className:"admonition-heading"},(0,r.kt)("h5",{parentName:"div"},(0,r.kt)("span",{parentName:"h5",className:"admonition-icon"},(0,r.kt)("svg",{parentName:"span",xmlns:"http://www.w3.org/2000/svg",width:"14",height:"16",viewBox:"0 0 14 16"},(0,r.kt)("path",{parentName:"svg",fillRule:"evenodd",d:"M7 2.3c3.14 0 5.7 2.56 5.7 5.7s-2.56 5.7-5.7 5.7A5.71 5.71 0 0 1 1.3 8c0-3.14 2.56-5.7 5.7-5.7zM7 1C3.14 1 0 4.14 0 8s3.14 7 7 7 7-3.14 7-7-3.14-7-7-7zm1 3H6v5h2V4zm0 6H6v2h2v-2z"}))),"important")),(0,r.kt)("div",{parentName:"div",className:"admonition-content"},(0,r.kt)("p",{parentName:"div"},"The Octopus package on Bioconda is built to a ",(0,r.kt)("a",{parentName:"p",href:"https://en.wikipedia.org/wiki/Haswell_(microarchitecture)"},"Haswell")," architecture. This means that it will only work on  Haswell (with AVX2) or newer machines. If you need another architecture then consider using ",(0,r.kt)("a",{parentName:"p",href:"https://docs.conda.io/projects/conda-build/en/latest/"},"conda-build"),"."))))}u.isMDXComponent=!0}}]);