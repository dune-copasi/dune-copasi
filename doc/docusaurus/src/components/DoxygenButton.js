
import Link from '@docusaurus/Link';
import React from 'react';
import clsx from 'clsx';
import useDocusaurusContext from '@docusaurus/useDocusaurusContext';

import { useActiveVersion } from '@docusaurus/plugin-content-docs/client'

export default function DoxygenButton() {
  const activeVersion = useActiveVersion();
  const doxygenPath = activeVersion.isLast ? "/api/"+activeVersion.label : activeVersion.path.replace("/docs/", "/api/")

  return (
    <div className="container">
      <h1 className="hero__title">{useDocusaurusContext().title}</h1>
      <p className="hero__subtitle">{useDocusaurusContext().tagline}</p>
      <Link
        className={clsx(
          "button button--outline button--secondary button--lg",
          "getStarted"
        )}
        to={`${doxygenPath}`}
      >
        Go to ({activeVersion.label}) Online Doxygen Documentation
      </Link>
    </div>
  );
}
