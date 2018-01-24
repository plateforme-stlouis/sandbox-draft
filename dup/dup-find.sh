#!/bin/bash --noprofile --norc

#set -e

DUP_HOME=$HOME/duplicates
DUP_DB=$DUP_HOME/all-indexed-dir.db

DUP_BS=1M
DUP_COUNT=10

# +0 means as many as possible
# (see GNU parallel, option -j)
DUP_NTHREADS=+0

DUP_GIT=0

DUP_TMP=$DUP_HOME/tmp-remove.that
DUP_LST=$DUP_HOME/lst-remove.that
DUP_HSH=$DUP_HOME/hsh-remove.that

function dup--rm-tmp () {
    [ -f $DUP_TMP ] && command rm -I $DUP_TMP
    [ -f $DUP_LST ] && command rm -I $DUP_LST
    [ -f $DUP_HSH ] && command rm -I $DUP_HSH
    return 0
}

function dup--init-maybe () {
    if [ ! -d $DUP_HOME ]
    then
        read -p "Do you want to initialize DUP? (yes/NO) " yn
        case $yn in
            y|yes|Yes|YES)
                dup-init
                return 0
                ;;
            n|no|No|NO)
                echo "Bye."
                return 1
                ;;
            *)
                echo "Wrong ! Try again..."
                dup--init-maybe
                ;;
        esac
    fi
    return 0
}

function dup-help () {
    echo "Usage: dup-cmd [options]"
    echo ""
    echo "o Init: dup-init"
    echo "o Index and store: dup-index dir-name [option]"
    echo "Only one directory is accepted per invocation."
    echo "  -p/--parallel: parallel hash (depend on GNU parallel)."
    echo "  --abspath: dir-name is an absolute path."
    echo "o List duplicates from store: dup-show"
    echo "o Which directories indexed: dup-ls-dirs"
    echo "o Compare MD5sum outputs: dup-compare out1 out2 out3 ... outN"
    return 0
}

function dup-init () {
    if [ ! -d $DUP_HOME ]
    then
        command mkdir -p $DUP_HOME/md5
        command mkdir -p $DUP_HOME/tps
        echo "Initialized empty DUP repository in $DUP_HOME"
    else
        if [ -d $DUP_HOME/.git ]
        then
            echo "It seems already initialized (with Git)."
        else
            echo "It seems already initialized (without Git)."
        fi
    fi
    if [ ! -d $DUP_HOME/.git ]
    then
        if [ ! $DUP_GIT -eq 0 ]
        then
            from=$(pwd)
            command cd $DUP_HOME
            git init
        fi
        touch $DUP_DB
        if [ ! $DUP_GIT -eq 0 ]
        then
            git add $DUP_DB
            git commit -m "init"
            command cd $from
        fi
    fi
    return 0
}

function dup-ls () {
    dup--init-maybe || return 0
    from=$(pwd)
    command cd $DUP_HOME/md5
    command find . -name "*.md5" -type f -print \
        | command sort | xargs -I {} basename {} .md5
    command cd $from
    return 0
}

function dup-ls-files () {
    dup--init-maybe || return 0
    if [ -d $DUP_HOME/.git/ ]
    then
        from=$(pwd)
        command cd $DUP_HOME
        git ls-files
        command cd $from
    else
        echo "Require Git."
        echo "Try: DUP_GIT=1 dup-gitify"
    fi
    return 0
}

function dup--db-add () {
    dup--init-maybe || return 0
    dup--rm-tmp
    if [ -f $DUP_DB ]
    then
        command cp $DUP_DB $DUP_TMP
    else
        command touch $DUP_TMP
    fi
    echo $1 >> $DUP_TMP
    command cat $DUP_TMP \
        | command sort   \
        | command uniq   \
                  > $DUP_DB
    dup--rm-tmp
    return 0
}

function dup-ls-dirs () {
    dup--init-maybe || return 0
    if [ -f $DUP_DB ]
    then
        cat $DUP_DB
    fi
    return 0
}

function dup-index () {
    dup--init-maybe || return 0
    local from=$(pwd)
    local DUP_PARALLEL=0
    local DUP_ABSPATH=0
    if [ $# -lt 1 ]
    then
        dup-help
        return 1
    fi
    input=$1
    shift
    while [ $# -gt 0 ]
    do
        case $1 in
            -p|--parallel)
                DUP_PARALLEL=1
                shift
                ;;
            --abspath)
                DUP_ABSPATH=1
                shift
                ;;
        esac
    done
    if [ $DUP_ABSPATH -eq 1 ]
    then
        abspath=$input
    else
        abspath=$from/$input
    fi
    if [ ! -d $abspath ]
    then
        echo "Error: $abspath does not exist."
        return 1
    fi

    echo "Indexing ${abspath}..."

    dup--rm-tmp
    if [ ! -d $DUP_HOME ]
    then
        dup-init
    fi
    # remove trailing slash /
    dir=${abspath%/}
    # replace / by -
    name="dir-${dir//\//-}"
    md5=md5/${name}.md5
    tps=tps/${name}.tps

    echo " o Hash (checksum)..."
    # cksum should hash faster, but it practically seems not.
    # but the bottleneck seems the IO: read the file!
    if [ $DUP_PARALLEL -eq 1 ]
    then
        echo "  ...Running in parallel..."
        find $dir ! -empty ! -type l -type f -print0 \
            | parallel --progress -j$DUP_NTHREADS -0 md5sum > $DUP_TMP
    else
        # find $dir ! -empty ! -type l -type f -print0 \
        #     | xargs -0 md5sum > $DUP_TMP
        find $dir ! -empty ! -type l -type f -print > $DUP_LST
        for file in $(cat $DUP_LST)
        do
            dd status=none if=$file bs=$DUP_BS count=$DUP_COUNT \
                | md5sum
        done | tr -s '-' ' ' > $DUP_HSH
        paste -d ' ' $DUP_HSH $DUP_LST > $DUP_TMP
    fi
    cat $DUP_TMP | sort > $DUP_HOME/$md5
    dup--rm-tmp

    echo " o Change time (stat)..."
    # This *UGLY* for-loop because quote in filename
    # and I did not find a way,
    # e.g., awk | xargs -0 which should fix the issue
    local is_error=0
    for file in $(awk '{print $2}' $DUP_HOME/$md5)
    do
        stat --printf="%Z %n\n" $file >> $DUP_HOME/$tps
        if [ ! $? -eq 0 ]
        then
            is_error=1
        fi
    done
    if [ $is_error -eq 1 ]
    then
        echo ""
        echo "Warning: filename containing quotes (' or \") or spaces."
        echo "Warning: you should rename it because it is a wrong habit."
    fi

    dup--db-add $dir

    if [ -d $DUP_HOME/.git/ ]
    then
        command cd $DUP_HOME
        git add $md5
        git add $tps
        git commit -am "$name"
        command cd $from
    fi

    return 0
}

function dup-log () {
    dup--init-maybe || return 0
    if [ -d $DUP_HOME/.git/ ]
    then
        from=$(pwd)
        command cd $DUP_HOME
        git --no-pager log --pretty="%Cblue%h %Creset%s"
        command cd $from
    else
        echo "Require Git."
        echo "Try: DUP_GIT=1 dup-gitify"
    fi
    return 0
}

function dup-reindex () {
    dup--init-maybe || return 0
    for dir in $(cat $DUP_DB)
    do
        dup-index $dir --abspath
    done
    return 0
}

function dup-gitify () {
    dup-init
    from=$(pwd)
    command cd $DUP_HOME
    for file in $(command ls -1 $DUP_HOME/md5/*.md5)
    do
        git add $file
        git commit -am "$file"
    done
    command cd $from
    return 0
}

function dup-show-raw (){
    dup--init-maybe || return 0
    dup-compare $(command ls -1 $DUP_HOME/md5/*.md5)
    return 0
}

function dup-show () {
    # Slow implementation !
    # But does the job ;-)
    local prev="a"
    for line in $(dup-show-raw | awk '{printf ("%s|%s\n", $1,$2)}')
    do
        id=$(echo $line | cut -d'|' -f1)
        file=$(echo $line | cut -d'|' -f2)
        [ ! "$id" == "$prev" ] && echo " "
        echo $id $file
        prev=$id
    done
}

function dup-compare () {
    for file in $*
    do
        cat $file
    done | sort | uniq -w32 -dD
    return 0
}

function dup--check-sanity () {
    dup--init-maybe || return 0

    # Not sure that parallel is faster
    # md5=$(ls -1 $DUP_HOME/md5/*.md5 | parallel cat | wc -l)
    # because the time of each cat

    local md5=0
    local tps=0
    md5=$(for file in $(command ls -1 $DUP_HOME/md5/*.md5)
          do
              cat $file
          done | wc -l)
    tps=$(for file in $(command ls -1 $DUP_HOME/tps/*.tps)
          do
              cat $file
          done | wc -l)
    if [ ! $md5 -eq $tps ]
    then
        echo "Error: $md5 MD5 files vs $tps change-time."
        echo "Error: Should coming from quotes (stat)."
        return 1
    else
        echo "Ok."
        return 0
    fi
}

function dup--check-sanity-par () {
    dup--init-maybe || return 0

    local md5=0
    local tps=0
    # md5=$(find $DUP_HOME/md5 -type f -print0 | parallel cat | wc -l)
    # tps=$(find $DUP_HOME/tps -type f -print0 | parallel cat | wc -l)
    md5=$(ls -1 $DUP_HOME/md5/*.md5 | parallel cat | wc -l)
    tps=$(ls -1 $DUP_HOME/tps/*.tps | parallel cat | wc -l)
    echo $md5 $tps
}


function dup-rm () {
    dup--init-maybe || return 0
    rm -i $DUP_HOME/md5/${1}.md5
    rm -i $DUP_HOME/tps/${1}.tps
    return 0
}

function dup-search-dir () {
    dup--init-maybe || return 0
    dup-ls-dirs | grep -e $1
}

function dup-search-file () {
    dup--init-maybe || return 0
    dup-ls | grep -e $1
}

function dup-version () {
    echo "It is too simple to have a version number."
    echo ""
    dup-help
}
