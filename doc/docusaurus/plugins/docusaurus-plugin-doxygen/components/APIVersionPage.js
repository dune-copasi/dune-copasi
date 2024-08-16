import React from "react";
import Layout from "@theme/Layout";

function APIVersionPage({version}) {
  return (
    <Layout title={`Doxygen Documentation - ${version}`}>
      <div style={{ height: "100vh" }}>
        <iframe
          src={`/doxygen/${version}/`}
          style={{ width: "100%", height: "100%", border: "none" }}
          title={`Doxygen Documentation ${version}`}
        />
      </div>
    </Layout>
  );
}

export default APIVersionPage;
