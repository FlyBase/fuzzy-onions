#!/bin/bash

SVM_DIR=$1
REGEX=$2

find $SVM_DIR -type f | while read FILE ; do
    BASENAME=${FILE##*/}
    case "$BASENAME" in
    FBrf*)
        echo -e "$BASENAME\t$FILE" >> $$.fulltext-by-fbrf
        ;;
    *)
        echo -e "$BASENAME\t$FILE" >> $$.fulltext-by-pmid
    esac
done

cat | while read FBRF PMID ; do
    FULLTEXT_PATH=$(cat $$.fulltext-by-fbrf | grep ^$FBRF | cut -f2)
    if [ -z "$FULLTEXT_PATH" ]; then
        FULLTEXT_PATH=$(cat $$.fulltext-by-pmid | grep ^$PMID | cut -f2)
    fi

    if [ -z "$FULLTEXT_PATH" ]; then
        echo -e "$FBRF\tFull-text not available"
    elif [ ! -f "$FULLTEXT_PATH" ]; then
        echo -e "$FBRF\tFull-text not available"
    else
        NMATCH=$(grep -i -c -E "$REGEX" $FULLTEXT_PATH)
        echo -e "$FBRF\t$NMATCH"
    fi
done

rm $$.fulltext-by-{fbrf,pmid}
