import React from 'react';
import clsx from 'clsx';
import styles from './HomepageFeatures.module.css';

const FeatureList = [
  {
    title: 'Accurate',
    // Svg: require('../../static/img/undraw_docusaurus_mountain.svg').default,
    description: (
      <>
        Octopus achieves class-leading accuracy in short-read small variant calling,
        including for germline and somatic mutations using a combination of Bayesian
        modelling and machine learning.
      </>
    ),
  },
  {
    title: 'Versatile',
    // Svg: require('../../static/img/undraw_docusaurus_tree.svg').default,
    description: (
      <>
        Designed to work with sequencing data generated from different 
        experimental setups, Octopus can call variants in germline, bulk tumour,
        metagenomic, and single-cell datasets. 
      </>
    ),
  },
  {
    title: 'Easy to use',
    // Svg: require('../../static/img/undraw_docusaurus_react.svg').default,
    description: (
      <>
        Common read pre-processing steps are handled internally. 
        A uniform interface simplifies workflows. Multithreading is built-in.
        Octopus is designed with the user in mind. 
      </>
    ),
  },
];

function Feature({Svg, title, description}) {
  return (
    <div className={clsx('col col--4')}>
      {/* <div className="text--center">
        <Svg className={styles.featureSvg} alt={title} />
      </div> */}
      <div className="text--center padding-horiz--md">
        <h3>{title}</h3>
        <p>{description}</p>
      </div>
    </div>
  );
}

export default function HomepageFeatures() {
  return (
    <section className={styles.features}>
      <div className="container">
        <div className="row">
          {FeatureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}
