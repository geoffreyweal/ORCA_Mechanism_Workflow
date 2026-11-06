#!/bin/bash
# Preview the documentation locally (http://localhost:8000). Needs: pip install zensical (Python 3.11 or newer)
cd "$(dirname "$0")"
zensical serve
