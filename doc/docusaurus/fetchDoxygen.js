const https = require('https');
const fs = require('fs');
const path = require('path');
const unzipper = require('unzipper');
const os = require('os');
const { promisify } = require('util');

const unlink = promisify(fs.unlink);

// Define the regex pattern to extract version
const versionRegex = /Version:\s*([\w.-]+)/;

const BASE_URL = 'https://gitlab.dune-project.org/api/v4/projects/819/packages/generic/dune-copasi/';

// Function to read text file and extract version
function readVersionFromFile(filePath) {
  return new Promise((resolve, reject) => {
    fs.readFile(filePath, 'utf8', (err, data) => {
      if (err) return reject(err);

      // Apply the regex pattern to find the version
      const match = data.match(versionRegex);
      if (match) {
        // Extract the version number
        resolve(match[1]);
      } else {
        // Handle the case where no version is found
        resolve(null);
      }
    });
  });
}

// Function to read JSON file
function readJsonFile(filePath) {
  return new Promise((resolve, reject) => {
    fs.readFile(filePath, 'utf8', (err, data) => {
      if (err) return reject(err);
      try {
        const json = JSON.parse(data);
        resolve(json);
      } catch (e) {
        reject(e);
      }
    });
  });
}


async function downloadAndUnzip(version, url) {
  // Create paths for temporary ZIP file and final extraction directory
  const tempZipPath = path.join(os.tmpdir(), `${version}.zip`);
  const versionDir = path.join(__dirname, 'static/doxygen', version);

  console.log(`Downloading ${url} to ${tempZipPath}`);

  // Download the Doxygen ZIP file to a temporary directory
  await new Promise((resolve, reject) => {
    const file = fs.createWriteStream(tempZipPath);
    https.get(url, (response) => {
      if (response.statusCode !== 200) {
        return reject(new Error(`Failed to download file. Status code: ${response.statusCode}`));
      }
      response.pipe(file);
      file.on('finish', () => {
        file.close(resolve);
      });
    }).on('error', (err) => {
      fs.unlink(tempZipPath, () => {}); // Ignore errors during unlink
      reject(err);
    });
  });

  // Check if the downloaded file contains an error message
  const fileContent = await fs.promises.readFile(tempZipPath, 'utf8').catch(() => ''); // Handle error if read fails
  if (fileContent.includes('404 Not Found')) {
    await unlink(tempZipPath); // Use unlink from promisify to handle errors
    console.error(`The file for version ${version} was not found on the server.`);
    return;
  }

  console.log(`Unzipping ${tempZipPath} to ${versionDir}`);

  // Create the extraction directory if it does not exist
  await fs.promises.mkdir(versionDir, { recursive: true });

  // Unzip the Doxygen ZIP file
  try {
    await fs.createReadStream(tempZipPath)
      .pipe(unzipper.Extract({ path: versionDir }))
      .promise();
  } catch (error) {
    console.error(`Error unzipping ${tempZipPath}:`, error);
    throw error;
  }

  // Clean up the downloaded ZIP file
  await unlink(tempZipPath); // Use unlink from promisify to handle errors
}

// Main function to process versions
async function main() {
  try {
    // Read previous versions from JSON file
    const previousVersions = await readJsonFile(path.join(__dirname, 'versions.json'));

    // Read the development version identifier from a text file
    const developmentVersion = await readVersionFromFile(path.join(__dirname, '../../dune.module'));

    // Construct URLs for previous versions
    const versions = previousVersions.reduce((acc, version) => {
      acc[version] = `${BASE_URL}${version.substring(1)}/dune-copasi-doxygen.zip`;
      return acc;
    }, {});

    // Add development version to the versions object
    versions['next'] = `${BASE_URL}${developmentVersion}/dune-copasi-doxygen.zip`;

    // Process all versions
    for (const [version, url] of Object.entries(versions)) {
      console.log(`Processing version ${version}...`);
      try {
        await downloadAndUnzip(version, url);
      } catch (error) {
        console.error(`Version ${version} could not be properly downloaded\n`, error);
      }
    }
    console.log('All versions processed and ready in static/doxygen');
  } catch (error) {
    console.error('Error processing versions:', error);
  } finally {
    // Ensure the process exits
    process.exit();
  }
}

// Run the main function
main().catch(error => {
  console.error('Uncaught Error:', error);
  process.exit(1);
});
