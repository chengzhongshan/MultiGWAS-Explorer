winget install --id GitHub.cli --source winget --scope user --silent --accept-package-agreements --accept-source-agreements
gh auth login --web --clipboard --git-protocol https --hostname github.com
gh auth setup-git --hostname github.com

GitHub CLI displayed a one-time code, copied it to the clipboard, and opened GitHub’s device-login page. You entered the code and authorized the login. No token needed to be pasted into the terminal or this chat. See the GitHub CLI login guide.
Because your Git is Cygwin Git, gh auth setup-git wrote a Windows path that Git could not execute. I corrected the GitHub credential helper like this:
$ghExe = (Get-Command gh).Source
$ghUnix = (& cygpath -u $ghExe).Trim()
$key = 'credential.https://github.com.helper'
git config --global --unset-all $key
git config --global --add $key ''
git config --global --add $key "!'$ghUnix' auth git-credential"
The empty helper entry stops Git from using the old stored credential; the next entry makes it ask GitHub CLI for the new one. For a future expired login, start with gh auth refresh --hostname github.com; you should not need to reinstall the CLI.
After testing and committing the code, I pushed it with:
cd 'C:\Users\zcheng\Downloads\AI4Coding\perlMCP4AOA_GWAS\MultiGWAS-Explorer'
git status -sb
git fetch origin
git diff --cached --check
git commit -m "Import combined meta-analysis GWAS directly in SAS and gnuplot workflows"
git push origin main
git status -sb