# Cloud pilot transfer

This optional technical pilot accompanies the refreshed 0.1.0 beta. Doug restored the pinned runtime and exercised all four tools on Posit Cloud Linux, with user-reported native agreement and resource measurements. The merged implementation passed 128 sequential MCP tests locally. Classroom adoption and scientific and Indonesian review remain pending. This ZIP contains the pinned scientific source package and no installed Mac libraries.

Upload `occJSDM-cloud-pilot.zip` into the Cloud project root using the Files pane. Posit Cloud may offer to unzip uploaded ZIP files. If it extracts automatically, retain the resulting `occJSDM-mcp` folder and skip manual extraction. Place that folder at `mcp/paper2agent`; do not overwrite an existing directory.

If the ZIP remains unextracted, run in the Cloud Terminal:

```sh
cd /cloud/project
unzip occJSDM-cloud-pilot.zip
mkdir -p mcp
test ! -e mcp/paper2agent && mv occJSDM-mcp mcp/paper2agent
```

Then install from the verified Python 3.12.11 available in this Cloud project:

```sh
cd /cloud/project/mcp/paper2agent
python3 -m venv occJSDM-env
./occJSDM-env/bin/python -m pip install -r src/requirements.txt
Rscript --vanilla r-runtime/restore-runtime.R "$PWD/cloud-r-runtime"
```

The R restore requires internet access, the native compiler/toolchain and enough RAM to compile occJSDM. It uses a new isolated library for the pinned package, so your existing ordinary installation remains separate. Save the installation output. If a step fails, inspect its error before continuing; `cloud-r-runtime` must be new or empty when starting restoration.

After restoration succeeds, merge the entry from `class/posit-assistant-settings.template.json` into your project's `.posit/assistant/settings.json`. Preserve existing entries. Trust this project workspace in Assistant, then inspect `/mcp`. The connected server should expose `validate_data`, `fit_model`, `diagnostics` and `summarise_fit`. Assistant adds a client namespace to the tool names. See [Posit's MCP configuration documentation](https://assistant.posit.co/docs/reference/mcp-servers/).

Do not paste an API key into this package. Model credentials stay in Assistant. See `USAGE.md` for tool requests, artifacts and scientific limits.
