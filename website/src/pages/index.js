import React from 'react';
import clsx from 'clsx';
import Layout from '@theme/Layout';
import Link from '@docusaurus/Link';
import useDocusaurusContext from '@docusaurus/useDocusaurusContext';
import styles from './index.module.css';
import HomepageFeatures from '../components/HomepageFeatures';
import useBaseUrl from '@docusaurus/useBaseUrl';

function HomepageHeader() {
  const {siteConfig} = useDocusaurusContext();
  return (
    <header className={clsx('hero hero--primary', styles.heroBanner)}>
      <div className="container">
        <h1 className="hero__title">
          <img src={useBaseUrl('/img/octopus_keytar.svg')} className={styles.heroLogo} alt="octopus" />
        </h1>
        <p className="hero__subtitle">{siteConfig.tagline}</p>
        <div className={styles.social}>
        <iframe 
          src="https://ghbtns.com/github-btn.html?user=luntergroup&amp;repo=octopus&amp;type=watch&amp;count=true" 
          height={30} 
          width={118} 
          frameBorder={0} 
          scrolling="0" 
          style={{ width: '118px', height: '30px' }}
          ></iframe>
        </div>
        <div className={styles.buttons}>
          <Link
            className="button button--secondary button--lg"
            to="/docs/installation">
          Install
          </Link>
          <Link
            className="button button--secondary button--lg"
            to="/docs/tutorials/germline">
          Use
          </Link>
          <Link
            className="button button--secondary button--lg"
            to="https://gitter.im/octopus/octopus">
          Ask
          </Link>
          <Link
            className="button button--secondary button--lg"
            to="/docs/publications">
          Cite
          </Link>
        </div>
      </div>
    </header>
  );
}

export default function Home() {
  const {siteConfig} = useDocusaurusContext();
  return (
    <Layout
      title={`${siteConfig.title}`}
      description="Haplotype-based variant calling">
      <HomepageHeader />
      <main>
        <HomepageFeatures />
      </main>
    </Layout>
  );
}
