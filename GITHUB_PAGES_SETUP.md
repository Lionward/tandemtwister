# GitHub Pages Setup Instructions

This repository now includes an `index.html` file with the Google site verification meta tag for indexing.

## How to Enable GitHub Pages

To make the Google site verification active, you need to enable GitHub Pages for this repository:

1. Go to your repository on GitHub: https://github.com/Lionward/tandemtwister
2. Click on **Settings** (in the repository menu)
3. In the left sidebar, click on **Pages**
4. Under "Build and deployment":
   - **Source**: Select "Deploy from a branch"
   - **Branch**: Select "main" (or your default branch) and "/ (root)"
   - Click **Save**

5. GitHub Pages will be deployed to: `https://lionward.github.io/tandemtwister/`

## After Enabling GitHub Pages

Once GitHub Pages is enabled:
1. The `index.html` file will be served at the GitHub Pages URL
2. Google will be able to verify your site using the meta tag in the `<head>` section
3. Visitors will be automatically redirected to the GitHub repository

## Verification

After enabling GitHub Pages, you can verify the setup by:
1. Visiting your GitHub Pages URL
2. Viewing the page source to confirm the meta tag is present
3. Using Google Search Console to complete the verification process

## Note

The current `index.html` includes:
1. Google site verification meta tag: `content="Fq4C_Fkmju44GXVeW2otbQB29bgU6nYfPrjctY2lEzg"`
2. Automatic redirect to the GitHub repository (3 second delay for accessibility)
