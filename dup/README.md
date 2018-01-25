
It is a dummy and really simple implementation to track duplicated files.

Just a bash wrapper around the command:
```bash
    find . ! -empty -type f -print0 \
        | xargs -0 md5              \
        | sort                      \
        | uniq -w32 -dD
```
nothing more. **And the bottleneck is the IO.**

The file MD5 indexes are stored under `$DUP_HOME` to browse duplicates.


And after the command `dup-update`, then the file `$DUP_DUP` contains the
duplicated files.

 - First pass: only `$DUP_BS*$DUP_COUNT`MB is hashed. (`dup-index`)
 - Second pass: the previous duplicated hashes are fully rehashed. (`dup-update`)


**Note:** is it faster than `dup-1` ? (`git co -b one dup-1`)


### Simple install
Fetch the `dup-find.sh` file and source it. Done :-)
```bash
    wget https://raw.githubusercontent.com/plateforme-stlouis/sandbox-draft/master/dup/dup-find.sh
    source dup-find.sh
    dup-version
```

### Example

Index your home directory:
```bash
    dup-init
    for dir in $(command ls -d1 $HOME/*) ; do dup-index $dir ; done
    dup-update # Should take some time

    dup-show
```

If GNU parallel is installed, then parallelize a bit some internal with the
option `-p`:
```bash
    for dir in $(command ls -1d $HOME/*)
    do
        dup-index $dir -p
    done

    dup-show
```
