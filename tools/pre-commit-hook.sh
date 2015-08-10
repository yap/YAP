ASTYLE=$(which astyle)
if [ $? -ne 0 ]; then
	echo "[commit rejected] astyle not installed. Unable to check source file format policy." >&2
	exit 1
fi

OPTIONS="--options=tools/astylerc --dry-run"

RETURN=0
git diff --cached --name-status --diff-filter=ACMR |
    {
	    # Command grouping to workaround subshell issues. When the while loop is
	    # finished, the subshell copy is discarded, and the original variable
	    # RETURN of the parent hasn't changed properly.
	    while read STATUS FILE; do
            # regex matching
	        if [[ "$FILE" =~ \.(C|cxx|h)$ ]]; then
		        formatted=`$ASTYLE $OPTIONS $FILE | sed -n '/^Unchanged/p'`
		        if [ -z "$formatted" ]; then
			        echo "[commit rejected] $FILE does not respect the agreed coding standards." >&2
			        RETURN=1
		        fi
	        fi
	    done

	    if [ $RETURN -eq 1 ]; then
		    echo "">&2
		    echo "Make sure to run 'make format' and review the changes *before* committing."  >&2
	    fi
	    exit $RETURN
    }