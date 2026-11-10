import sys
from data.fetch_chembl import fetch_dpp4_data

if __name__ == "__main__":
    try:
        fetch_dpp4_data()
    except RuntimeError as err:
        print(str(err), file=sys.stderr)
        sys.exit(1)
