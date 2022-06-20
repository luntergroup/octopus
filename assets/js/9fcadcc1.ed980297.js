"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[81],{3905:function(e,t,n){n.d(t,{Zo:function(){return m},kt:function(){return c}});var r=n(7294);function a(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function i(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);t&&(r=r.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,r)}return n}function o(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?i(Object(n),!0).forEach((function(t){a(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):i(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function l(e,t){if(null==e)return{};var n,r,a=function(e,t){if(null==e)return{};var n,r,a={},i=Object.keys(e);for(r=0;r<i.length;r++)n=i[r],t.indexOf(n)>=0||(a[n]=e[n]);return a}(e,t);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);for(r=0;r<i.length;r++)n=i[r],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(a[n]=e[n])}return a}var s=r.createContext({}),p=function(e){var t=r.useContext(s),n=t;return e&&(n="function"==typeof e?e(t):o(o({},t),e)),n},m=function(e){var t=p(e.components);return r.createElement(s.Provider,{value:t},e.children)},f={inlineCode:"code",wrapper:function(e){var t=e.children;return r.createElement(r.Fragment,{},t)}},u=r.forwardRef((function(e,t){var n=e.components,a=e.mdxType,i=e.originalType,s=e.parentName,m=l(e,["components","mdxType","originalType","parentName"]),u=p(n),c=a,d=u["".concat(s,".").concat(c)]||u[c]||f[c]||i;return n?r.createElement(d,o(o({ref:t},m),{},{components:n})):r.createElement(d,o({ref:t},m))}));function c(e,t){var n=arguments,a=t&&t.mdxType;if("string"==typeof e||a){var i=n.length,o=new Array(i);o[0]=u;var l={};for(var s in t)hasOwnProperty.call(t,s)&&(l[s]=t[s]);l.originalType=e,l.mdxType="string"==typeof e?e:a,o[1]=l;for(var p=2;p<i;p++)o[p]=n[p];return r.createElement.apply(null,o)}return r.createElement.apply(null,n)}u.displayName="MDXCreateElement"},4775:function(e,t,n){n.r(t),n.d(t,{default:function(){return f},frontMatter:function(){return l},metadata:function(){return s},toc:function(){return p}});var r=n(7462),a=n(3366),i=(n(7294),n(3905)),o=["components"],l={id:"forest",title:"Random Forest"},s={unversionedId:"guides/filtering/forest",id:"guides/filtering/forest",isDocsHomePage:!1,title:"Random Forest",description:"Octopus provides a powerful way to classify variant calls with random forests using the Ranger random forest library. Pre-trained random forests are available on Google Cloud. Currently there are forests for germline and somatic variants. You can easily obtain the forests using the Python install script:",source:"@site/docs/guides/filtering/forest.md",sourceDirName:"guides/filtering",slug:"/guides/filtering/forest",permalink:"/octopus/docs/guides/filtering/forest",editUrl:"https://github.com/${organizationName}/${projectName}/edit/${branch}/website/docs/guides/filtering/forest.md",version:"current",frontMatter:{id:"forest",title:"Random Forest"},sidebar:"docs",previous:{title:"Introduction",permalink:"/octopus/docs/guides/filtering/annotations"},next:{title:"Hard Thresholds",permalink:"/octopus/docs/guides/filtering/thresholds"}},p=[{value:"Using random forest filtering",id:"using-random-forest-filtering",children:[{value:"Germline variant random forest filtering",id:"germline-variant-random-forest-filtering",children:[]},{value:"Somatic variant random forest filtering",id:"somatic-variant-random-forest-filtering",children:[]},{value:"Re-filtering an Octopus VCF with random forests",id:"re-filtering-an-octopus-vcf-with-random-forests",children:[]}]},{value:"Training random forests",id:"training-random-forests",children:[]}],m={toc:p};function f(e){var t=e.components,n=(0,a.Z)(e,o);return(0,i.kt)("wrapper",(0,r.Z)({},m,n,{components:t,mdxType:"MDXLayout"}),(0,i.kt)("p",null,"Octopus provides a powerful way to classify variant calls with random forests using the ",(0,i.kt)("a",{parentName:"p",href:"https://github.com/imbs-hl/ranger"},"Ranger")," random forest library. Pre-trained random forests are available on ",(0,i.kt)("a",{parentName:"p",href:"https://console.cloud.google.com/storage/browser/luntergroup/octopus/forests/?project=parabolic-eon-208710"},"Google Cloud"),". Currently there are forests for germline and somatic variants. You can easily obtain the forests using the Python install script:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-shell"},"$ ./scripts/install.py --forests\n")),(0,i.kt)("p",null,"will download the forests into ",(0,i.kt)("inlineCode",{parentName:"p"},"octopus/resources/forests/"),". The provided forests are suitable for calling typical diploid germline and cancer data. They ",(0,i.kt)("em",{parentName:"p"},"may")," work well for other types of samples but this is untested. "),(0,i.kt)("h2",{id:"using-random-forest-filtering"},"Using random forest filtering"),(0,i.kt)("h3",{id:"germline-variant-random-forest-filtering"},"Germline variant random forest filtering"),(0,i.kt)("p",null,"To filter germline variants using the germline random forest, just specify the path to the forest in the ",(0,i.kt)("inlineCode",{parentName:"p"},"--forest-file")," option:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-shell"},"$ octopus -R hs37d5.fa -I NA12878.bam \\\n    --forest resources/forests/germline.forest\n")),(0,i.kt)("p",null,"All calls will be annotated with the ",(0,i.kt)("inlineCode",{parentName:"p"},"FORMAT")," field ",(0,i.kt)("inlineCode",{parentName:"p"},"RFGQ"),", which is the phred-scaled probability of the ",(0,i.kt)("inlineCode",{parentName:"p"},"GT")," call being correct according to the forest model. Each record will also be annotated with an ",(0,i.kt)("inlineCode",{parentName:"p"},"INFO")," field, ",(0,i.kt)("inlineCode",{parentName:"p"},"RFGQ_ALL"),", which is the phred-scalled probability that all ",(0,i.kt)("inlineCode",{parentName:"p"},"GT")," calls in the record are correct, this probability is used to filter calls (see ",(0,i.kt)("a",{parentName:"p",href:"/octopus/docs/cli#--min-forest-quality"},(0,i.kt)("inlineCode",{parentName:"a"},"--min-forest-quality")),") with the ",(0,i.kt)("inlineCode",{parentName:"p"},"RF")," filter."),(0,i.kt)("h3",{id:"somatic-variant-random-forest-filtering"},"Somatic variant random forest filtering"),(0,i.kt)("p",null,"When calling germline and somatic variants, you need to provide both random forests:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-shell"},"$ octopus -R hs37d5.fa -I normal.bam tumour.bam -N NORMAL \\\n    --forest resources/forests/germline.forest \\\n    --somatic-forest resources/forests/somatic.forest\n")),(0,i.kt)("p",null,"However, if you are calling somatic variants only (i.e. using the ",(0,i.kt)("inlineCode",{parentName:"p"},"--somatics-only")," flag), then you just need to provide the somatic forest:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-shell"},"$ octopus -R hs37d5.fa -I normal.bam tumour.bam -N NORMAL \\\n    --somatics-only \\\n    --somatic-forest resources/forests/somatic.forest\n")),(0,i.kt)("h3",{id:"re-filtering-an-octopus-vcf-with-random-forests"},"Re-filtering an Octopus VCF with random forests"),(0,i.kt)("p",null,"Random forest filtering can be applied to a VCF file produced by Octopus without recalling with the ",(0,i.kt)("inlineCode",{parentName:"p"},"--filter-vcf")," option:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-shell"},"$ octopus -R hs37d5.fa -I NA12878.bam \\\n    --filter-vcf unfiltered.vcf.gz \\\n    --forest resources/forests/germline.forest \\\n    -o filtered.vcf.gz\n")),(0,i.kt)("p",null,"Note this will overwrite any existing filtering information. Note, the VCF file cannot be a ",(0,i.kt)("inlineCode",{parentName:"p"},"--legacy")," VCF file, and should contain all variant calls that Octopus would make by default (i.e. do not use a VCF produced with the ",(0,i.kt)("inlineCode",{parentName:"p"},"--somatics-only")," option)."),(0,i.kt)("h2",{id:"training-random-forests"},"Training random forests"),(0,i.kt)("p",null,"The supplied random forest should perform well for the majority of users, however, it is possible that performance can be improved by training a new forest model on your own data."),(0,i.kt)("p",null,"In addition to Octopus requirements, training new forest models requires the following:"),(0,i.kt)("ol",null,(0,i.kt)("li",{parentName:"ol"},"A truth set of variants and high-confidence regions (e.g., GIAB or SynDip)."),(0,i.kt)("li",{parentName:"ol"},(0,i.kt)("a",{parentName:"li",href:"https://snakemake.readthedocs.io/en/stable/"},"Snakemake"),"."),(0,i.kt)("li",{parentName:"ol"},(0,i.kt)("a",{parentName:"li",href:"https://www.realtimegenomics.com/products/rtg-tools"},"RTG Tools"),".")),(0,i.kt)("p",null,"Forest models can - and ideally should - be trained on multiple examples runs (i.e., a run of Octopus). Each example can itself be calls for a single sample, or joint calls for multiple samples. For joint calls, each sample is used to generate training data independently, and a subset of samples can be selected."),(0,i.kt)("p",null,"Training new forests is done using the bundled Snakemake script ",(0,i.kt)("a",{parentName:"p",href:"https://github.com/luntergroup/octopus/blob/develop/scripts/forest.smk"},"forest.smk"),". Basic usage is like:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-shell"},"$ snakemake \\\n    --snakefile forest.smk \\\n    --configfile config.yaml \\\n    --cores 20\n")),(0,i.kt)("p",null,"You can of course use any ",(0,i.kt)("a",{parentName:"p",href:"https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html"},"Snakemake option"),". Note that the RTG Tools binary ",(0,i.kt)("inlineCode",{parentName:"p"},"rtg")," must be in your path."),(0,i.kt)("p",null,"The configuration file (YAML or JSON format) provided to Snakemake specifies the training data and setup. The format of the config file is:"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"truths")," is a dictionary of truth sets, the key being some unique label (e.g. ",(0,i.kt)("inlineCode",{parentName:"li"},"SAMPLE1.truth"),"). The label does ",(0,i.kt)("strong",{parentName:"li"},"not")," need to correspond to the sample name. Each truth set contains a sub dictionary with ",(0,i.kt)("inlineCode",{parentName:"li"},"vcf")," and ",(0,i.kt)("inlineCode",{parentName:"li"},"bed")," value pairs corresponding to the truth variants and high confidence regions, respectively."),(0,i.kt)("li",{parentName:"ul"},"examples is a list of examples to use. Each example is dictionary with the following fields:",(0,i.kt)("ul",{parentName:"li"},(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"name"),": the filename prefix to use for this example. If this is not specified then the name is automatically set by concatenating input BAM names. However, this can result in long filenames that may cause OS errors if many inputs are used."),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"reference"),": the reference fasta file. ","[",(0,i.kt)("strong",{parentName:"li"},"required"),"]"),(0,i.kt)("li",{parentName:"ul"},"'reads': A list of BAM files to use. The list may be omitted if a single BAM is used. ","[default none]"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"truth"),": Either the label of the truth set for single sample calling, or a dictionary of sample-truth pairs, where the key is the sample name (in the BAM file), and the value is one of the keys in ",(0,i.kt)("inlineCode",{parentName:"li"},"truths"),". ","[",(0,i.kt)("strong",{parentName:"li"},"required"),"]"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"regions"),": A bed file to use for calling (i.e. sets option ",(0,i.kt)("inlineCode",{parentName:"li"},"--regions-file"),"). ","[default none]"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"options"),": Additional options to pass to the octopus command."),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"tp_fraction"),": The fraction of true positive calls to use for training. ","[default 1]"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"fp_fraction"),": The fraction of false positive calls to use for training. ","[default 1]"),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"threads"),": How many threads to use for analysis in this run."))),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"training")," is a dictionary of forest training setting:",(0,i.kt)("ul",{parentName:"li"},(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"training_fraction"),": the fraction of calls to use for training."),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"hyperparameters"),":the random forest hyperparameters to use for training:",(0,i.kt)("ul",{parentName:"li"},(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"trees"),": the number of trees to use (ranger option ",(0,i.kt)("inlineCode",{parentName:"li"},"--ntrees"),")."),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"min_node_size"),": minimum size of a node in the tree (ranger option ",(0,i.kt)("inlineCode",{parentName:"li"},"--targetpartitionsize"),")."),(0,i.kt)("li",{parentName:"ul"},(0,i.kt)("inlineCode",{parentName:"li"},"maxdepth"),": maximum depth of each tree (ranger option ",(0,i.kt)("inlineCode",{parentName:"li"},"--maxdepth"),").")))))),(0,i.kt)("p",null,"A minimal example is:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-yaml"},"examples:\n    -\n        reference: /path/to/reference.fa\n        reads: /path/to/reads.bam\n        truth: GIAB//GRCh38//HG001\n")),(0,i.kt)("p",null,"Here, ",(0,i.kt)("inlineCode",{parentName:"p"},"truth")," is set to a special label ",(0,i.kt)("inlineCode",{parentName:"p"},"GIAB//GRCh38//HG001"),". The script accepts special labels like this for any of the GIAB truth sets - it will download all the required truth data. A more detailed example is: "),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-yaml"},"truths:\n    SAMPLE1.truth:\n        vcf: /path/to/truth/variants/sample1.vcf.gz\n        bed: /path/to/truth/variants/sample1.bed\n    SAMPLE2.truth:\n        vcf: /path/to/truth/variants/sample2.vcf.gz\n        bed: /path/to/truth/variants/sample2.bed\nexamples:\n    -\n        reference: /path/to/reference.fa\n        reads: /path/to/reads1.bam\n        truth: SAMPLE1.truth\n        tp_fraction: 0.5\n        fp_fraction: 1\n        threads: 10\n    -\n        name: joint\n        reference: /path/to/reference.fa\n        reads:\n            - /path/to/reads1.bam\n            - /path/to/reads2.bam\n        regions: /path/to/calling/regions.bed\n        truth:\n            SAMPLE1: SAMPLE1.truth\n            SAMPLE2: SAMPLE2.truth\n        options: --config /path/to/octopus/octopus.config\n        threads: 20\ntraining:\n    training_fraction: 0.25\n    hyperparameters:\n        -\n            trees: 200\n            min_node_size: 20\n")),(0,i.kt)("p",null,"The default behaviour of the script is to use all calls for training, however if the option ",(0,i.kt)("inlineCode",{parentName:"p"},"--kind=somatic")," is used then only variants called ",(0,i.kt)("inlineCode",{parentName:"p"},"SOMATIC")," are used for training. This can be used to generate somatic forests."))}f.isMDXComponent=!0}}]);