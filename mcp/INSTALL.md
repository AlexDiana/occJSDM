# Installing the occJSDM assistant (MCP) server

This page gets the optional occJSDM MCP server running under Posit Assistant, either on your own Mac or in a Posit Cloud project. You do not need it to use the package. For what the four tools do, how artifacts are written and what has been tested, see the [server usage guide](paper2agent/USAGE.md).

Cloning this repository does not install the server. The server needs a Python environment and its own R runtime, built from the occJSDM 0.3.0 source tarball, and Git carries neither. Everything is in one download, `occJSDM-cloud-pilot-v0.3.0.zip`, attached to the [v0.3.0-beta release](https://github.com/AlexDiana/occJSDM/releases/tag/v0.3.0-beta).

## On your own Mac

You install the server once, into its own folder outside any Git repository. Any number of projects, such as your occJSDM clone or a workshop folder, can then use that one installation. Each project only needs an entry in its own `.posit/assistant/settings.json` that points at it. Installing never touches your projects' tracked files.

These steps follow the restore recorded on macOS 27 (arm64) with R 4.6.1. The restore took about 7 minutes on an Apple silicon Mac with a fast connection, mostly installing 60 R packages, so allow 10 to 30 minutes. Keep using the same Terminal tab from step 2 onwards.

1.  Check the prerequisites. You need R 4.6.1, Python 3.11 or later, the Xcode command line tools, a Fortran compiler for R, and Homebrew's `gettext`.

    ``` sh
    Rscript --version
    python3 --version
    ls /opt/homebrew/opt/gettext/include
    ```

2.  Download the zip and unpack it. The unpacked folder is called `occJSDM-mcp`. The commands put it in your home folder, so replace `~` after `-d` with another folder, such as `~/src`, if you prefer. The folder must not exist already.

    ``` sh
    cd ~/Downloads
    curl -fL -O https://github.com/AlexDiana/occJSDM/releases/download/v0.3.0-beta/occJSDM-cloud-pilot-v0.3.0.zip
    unzip -q occJSDM-cloud-pilot-v0.3.0.zip -d ~
    cd ~/occJSDM-mcp
    ```

3.  Create the Python environment.

    ``` sh
    python3 -m venv occJSDM-env && ./occJSDM-env/bin/python -m pip install -r src/requirements.txt
    ```

4.  Restore the R runtime into a new folder. The two exports are the macOS workarounds. `MACOSX_DEPLOYMENT_TARGET` stops the Fortran compiler passing clang an invalid minimum macOS version, and the throwaway Makevars file lets the `rpart` build find Homebrew's gettext headers. That file replaces your own `~/.R/Makevars` for this Terminal session only. The restore checks the occJSDM source tarball against its pinned hash, so a damaged download is rejected. The target folder must not exist yet, so after a failed attempt delete `local-r-runtime` before retrying.

    ``` sh
    export MACOSX_DEPLOYMENT_TARGET=15.0
    export R_MAKEVARS_USER="$(mktemp)"
    echo 'CPPFLAGS += -I/opt/homebrew/opt/gettext/include' > "$R_MAKEVARS_USER"
    Rscript --vanilla r-runtime/restore-runtime.R "$PWD/local-r-runtime" 2>&1 | tee restore.log
    ```

    The last line should say that the restore and its isolated origin and version checks passed.

5.  Print the settings entry for a project. Set `PROJECT` to the folder you open in RStudio, then run the block from `~/occJSDM-mcp`. The server writes its results to an `artifacts` folder inside `PROJECT`, and it can only read data files under `PROJECT`.

    ``` sh
    PROJECT="$HOME/src/occJSDM"
    cat <<EOF
    {
      "mcpServers": {
        "occJSDM": {
          "type": "local",
          "command": ["$PWD/occJSDM-env/bin/python", "$PWD/src/occJSDM_mcp.py"],
          "environment": {
            "P2A_PROJECT_ROOT": "$PROJECT",
            "P2A_R_PROJECT": "$PWD/local-r-runtime",
            "P2A_RSCRIPT": "$(which Rscript)",
            "P2A_TIMEOUT_SECONDS": "600"
          },
          "timeout": 30000
        }
      }
    }
    EOF
    ```

6.  Merge the printed `occJSDM` entry into `PROJECT/.posit/assistant/settings.json`, creating the file if it does not exist and keeping any keys already there. Both the occJSDM and occJSDM-workshop repositories already tell Git to ignore that file and the `artifacts` folder. In any other project, check that Git ignores them before you commit. To use the server from a second project, change `PROJECT` and repeat steps 5 and 6. Do not install again.

7.  Restart the Posit Assistant session, trust the workspace when asked, and enter your model provider's key in Assistant's own credential settings, never in the project. Then follow the check at the end of this page.

## On Posit Cloud

Posit Cloud is Linux, so nothing needs the Mac workarounds. The restore took about 13 minutes on Cloud.

1.  Upload `occJSDM-cloud-pilot-v0.3.0.zip` from the release page to the root of the Cloud project with the Files pane. Cloud unpacks it into a folder called `occJSDM-mcp`. Delete the zip afterwards.

2.  In the Terminal tab, move the folder into place, then create the Python environment and restore the R runtime.

    ``` sh
    cd /cloud/project && mkdir -p mcp && mv occJSDM-mcp mcp/paper2agent
    ```

    ``` sh
    cd /cloud/project/mcp/paper2agent && python3 -m venv occJSDM-env && ./occJSDM-env/bin/python -m pip install -r src/requirements.txt
    ```

    ``` sh
    cd /cloud/project/mcp/paper2agent && Rscript --vanilla r-runtime/restore-runtime.R "$PWD/cloud-r-runtime" 2>&1 | tee restore.log
    ```

3.  Register the server. If the project has no `.posit/assistant/settings.json` yet, copy the template. Otherwise add its `occJSDM` entry to the existing `mcpServers` section.

    ``` sh
    mkdir -p /cloud/project/.posit/assistant && cp /cloud/project/mcp/paper2agent/class/posit-assistant-settings.template.json /cloud/project/.posit/assistant/settings.json
    ```

4.  Trust the workspace in Posit Assistant, choose a model, and enter the provider's key in Assistant's credential settings. Then follow the check below.

## Checking that it connected

Send `/mcp` in Posit Assistant. The Session information panel should have an MCP Servers section listing occJSDM as Connected. The panel does not list tools, so ask the model which occJSDM tools it has. It should name exactly `validate_data`, `fit_model`, `diagnostics` and `summarise_fit`, shown as `mcp__occJSDM__validate_data` and so on.

If the panel has no MCP Servers section at all, Assistant did not load any server:

- The `occJSDM` entry is missing from the `.posit/assistant/settings.json` of the project you actually opened.
- The workspace is not trusted. Assistant silently omits project-defined servers from an untrusted workspace.
- The settings file is not valid JSON, for example a missing comma after merging.

If occJSDM is listed but not Connected, a path in the entry is wrong or the runtime is missing. Run the server by hand to see the error. A healthy server prints nothing and waits for input, so stop it with Control-C. On Cloud, use the `/cloud/project/mcp/paper2agent` folder and `cloud-r-runtime`.

``` sh
cd ~/occJSDM-mcp
P2A_R_PROJECT="$PWD/local-r-runtime" ./occJSDM-env/bin/python src/occJSDM_mcp.py
```
