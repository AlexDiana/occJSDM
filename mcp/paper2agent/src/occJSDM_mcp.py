"""Stdio entry point for the independently verified occJSDM pilot tools."""
from fastmcp import FastMCP
from tools.quickstart import quickstart_mcp

mcp = FastMCP("occJSDM")
mcp.mount(quickstart_mcp)

if __name__ == "__main__":
    mcp.run(transport="stdio", show_banner=False)
