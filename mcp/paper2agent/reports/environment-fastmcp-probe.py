import asyncio
import json
from fastmcp import Client, FastMCP

server = FastMCP("Paper2MCP environment probe")

@server.tool()
def environment_echo(value: str) -> dict:
    return {"value": value}

async def main():
    async with Client(server) as client:
        tools = await client.list_tools()
        result = await client.call_tool("environment_echo", {"value": "isolated-runtime"})
        assert [tool.name for tool in tools] == ["environment_echo"]
        assert result.data == {"value": "isolated-runtime"}
        print(json.dumps({"inventory": [tool.name for tool in tools], "result_data": result.data, "status": "passed"}))

asyncio.run(main())
