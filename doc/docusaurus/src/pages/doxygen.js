import React from 'react';
import { useLocation } from 'react-router-dom';
import Layout from '@theme/Layout';

function DoxygenPage() {
  const { search } = useLocation();
  const params = new URLSearchParams(search);
  const version = params.get('page') || 'next';

  // Construct the URL for the Doxygen documentation based on the version
  const doxygenPath = `/doxygen/${version}/index.html`;

  return (
    <Layout title={`Doxygen Documentation - ${version}`}>
      <div style={{ height: '100vh' }}>
        <iframe
          src={doxygenPath}
          style={{ width: '100%', height: '100%', border: 'none' }}
          title={`Doxygen Documentation ${version}`}
        />
      </div>
    </Layout>
  );
}

export default DoxygenPage;
