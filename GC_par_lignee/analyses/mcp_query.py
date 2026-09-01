#!/usr/bin/env python3
"""
Objet       : interroger un serveur MCP streamable-HTTP (ici TBannotator) depuis
              un script, sans passer par l'agent. Le resultat va directement dans
              un fichier au lieu de transiter par le contexte : une extraction de
              plusieurs centaines de lignes de metadonnees ne coute alors rien.
              Protocole : POST initialize, notification initialized, puis
              tools/call ; la reponse arrive en SSE (ligne "data:").
Entrees     : requete SQL sur stdin
Sorties     : CSV sur stdout
Usage       : echo "SELECT ..." | python3 mcp_query.py > sortie.csv
Reutilisable: oui -- tout serveur MCP HTTP du poste, changer URL et TOOL
Projet      : GC_par_lignee
Date        : 2026-08-29
"""
import ast
import json
import sys
import urllib.request

URL = "https://tblearn.tbannotator.ideev.universite-paris-saclay.fr/mcp"
TOOL = "tool_query_postgres"
HEADERS = {"Content-Type": "application/json",
           "Accept": "application/json, text/event-stream"}


def post(payload, sid=None, timeout=600):
    h = dict(HEADERS)
    if sid:
        h["Mcp-Session-Id"] = sid
    req = urllib.request.Request(URL, json.dumps(payload).encode(), h)
    r = urllib.request.urlopen(req, timeout=timeout)
    return r.headers.get("Mcp-Session-Id"), r.read().decode()


def parse_sse(body):
    for line in body.splitlines():
        if line.startswith("data:"):
            return json.loads(line[5:].strip())
    return json.loads(body)


def query(sql, max_rows=5000, timeout=300):
    sid, _ = post({"jsonrpc": "2.0", "id": 1, "method": "initialize",
                   "params": {"protocolVersion": "2025-06-18",
                              "capabilities": {},
                              "clientInfo": {"name": "mcp_query", "version": "1"}}})
    post({"jsonrpc": "2.0", "method": "notifications/initialized"}, sid)
    _, body = post({"jsonrpc": "2.0", "id": 2, "method": "tools/call",
                    "params": {"name": TOOL,
                               "arguments": {"query": sql, "max_rows": max_rows,
                                             "timeout_seconds": timeout}}}, sid)
    txt = parse_sse(body)["result"]["content"][0]["text"]
    d = ast.literal_eval(txt) if txt.lstrip().startswith("{'") else json.loads(txt)
    if isinstance(d, str):
        d = ast.literal_eval(d)
    if not d.get("success"):
        raise RuntimeError(json.dumps(d)[:2000])
    if d["data"].get("truncated"):
        print("# ATTENTION : resultat tronque par max_rows", file=sys.stderr)
    return d["data"]["csv"]


if __name__ == "__main__":
    print(query(sys.stdin.read()), end="")
