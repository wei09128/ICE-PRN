# Publish after choosing the destination repository

Extract this ZIP as its own folder under `/mnt/c/Wei/ICE-PRN`. Do not overlay it onto an existing repository or the active research workspace. Verify and rebuild first:

```bash
cd /mnt/c/Wei/ICE-PRN/ICE_PRN_GITHUB_RELEASE_v1
python verify_release.py
python figure_runner/build_all.py
```

For a new, empty GitHub repository, use the real repository URL below. Authentication and your Git author identity must already be configured. No force push is needed.

```bash
git init -b main
git add README.md PUBLISH.md .gitignore requirements.txt verify_release.py RELEASE_HASHES.json VALIDATION.json docs paper figure_runner toc analysis_packages
git commit -m "Add ICE-PRN manuscript reproduction snapshot"
git remote add origin YOUR_GITHUB_REPOSITORY_URL
git push -u origin main
```

For an existing repository, supply its URL first so the update can be prepared against its current branch while preserving its files and history. Do not run the new-repository commands in that checkout.

Before calling this an open-source release, choose the code license. Before citing a public release in the paper, replace repository/DOI placeholders with the actual published URLs. The cover letter is intentionally excluded.
