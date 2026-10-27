# Rebuilding the repository without the tracked sequence data

**Decision (26 October 2026):** commit history is not worth preserving at this stage, so the
repository is rebuilt from the current working tree rather than rewritten. This removes the
need for `git filter-repo`, for `git fsck`, and for resolving the OneDrive placeholder
problem that blocks history traversal - none of those apply when the old history is simply
discarded.

**The sequence data never moves.** Every step operates on `.git`, not on the working tree,
so `Data/00-AAC-Sequences/CleanBactSeqs/` and
`Data/01-AAC-AssignSeqVariants/DADAFilteredSeqs/` stay where they are and the pipelines
keep running. No backup-and-restore is required.

**All commands below are Windows cmd.** PowerShell and bash equivalents are in the table at
the end.

---

## Starting state, verified 27 October 2026

Re-checked immediately before this run, read-only, from the repo root:

| | |
|---|---|
| `.git` | directory, 3,296 files, 7,504 MB, `HEAD` = `3372f483` on `main`, `origin` intact |
| `_git-old` | absent - yesterday's rollback moved it back into place |
| `_git-dead` | **present**, 4 files, 0 MB, `ReadOnly` - the husk from yesterday |
| `Data/00-AAC-Sequences/_mistrimmed_superseded` | 44 files, 112 MB, confirmed ignored |
| `Data/00-AAC-Sequences/CleanBactSeqs` | 322 files |

**`_git-dead` already exists, so do not use that name again.** The Recovery section below
tells you to `move .git _git-dead`; with the husk still in place `move` will fail or nest the
new directory inside the old one. Delete the husk first - `rm -rf _git-dead` in MobaXterm, or
`rmdir /s /q _git-dead` from your elevated cmd - or pick a fresh name such as `_git-dead2`.

**The optional trim has already been applied** (27 October). You are about to re-knit steps 01-05, which regenerates
`01-AAC-AssignSeqVariants_dadaF.ASV.rds` (37.6 MB) and `FastTree.log` (17.2 MB). Committing
them now and again after the re-run writes both versions into history permanently, which works
against the whole point of the exercise. Those three entries are now in `.gitignore`, which takes the first commit from 170 files /
134 MB to **158 / 79 MB** and LFS from 57 MB to **19 MB**. Re-verified by simulation.

---

## Outcome, 27 October 2026

**The orphan-commit route is the one that worked.** `move .git _git-old` never succeeded, so
steps 2 and 3 of the original route were never completed; everything below the heading
*Alternative route: orphan commit* is the procedure that was actually executed.

| | Before | After |
|---|---|---|
| `.git` total | 7,504 MB | **368 MB** |
| `.git/objects` | ~2,870 MB of packs plus loose objects | 349 MB |
| `.git/lfs` | 4,222 MB | **18 MB** |
| GitHub LFS quota | ~1.5 GB | ~19 MB |
| Tracked files | 618, incl. 479 `.fastq.gz` | 158, no FASTQ |
| LFS objects | 1,768 | 19 |
| History | 5 commits | `46e1692`, then `a8a8458` for this runbook |

`Data/00-AAC-Sequences/CleanBactSeqs` still holds its 322 files and
`_mistrimmed_superseded` its 44; no sequence data moved at any point.

### What each cleanup command actually did

* `git gc --prune=now` **works** on the object store - it took `.git/objects` from roughly
  2.9 GB to 349 MB - but its final cleanup sweep enters an unbounded per-directory
  *"Deletion of directory failed. should I try again (y/n)?"* loop. Answer `n` a few times to
  confirm the pattern, then **Ctrl+C**. This is safe: gc writes the new pack before deleting
  anything, so an interrupted prune leaves more objects than necessary, never fewer.
* An interrupted gc leaves a stale `.git/gc.pid`. Delete it (`del .git\gc.pid`) or later gc
  runs refuse with *gc is already running on machine …*.
* **The remaining bulk was never git objects.** At 4,571 MB total, 4,222 MB of it was the
  `.git/lfs` cache, which `git gc` does not touch. `git lfs prune` is the command that clears
  it, and because it deletes *files* rather than directories it is not blocked by the
  delete-child DENY ACE - it ran cleanly and reclaimed 4.2 GB in one pass. Run it before
  concluding that the space cannot be recovered.

### Two gotchas worth remembering

* `git rm -r --cached .` fails on an orphan branch with *"files have staged content different
  from both the file and the HEAD"*, listing any file edited since the last commit. Use
  `git rm -r --cached -f .`. With `--cached`, `-f` cannot touch the working tree; it only
  overrides a check that protects content existing solely in the index, and here the content
  is in the working tree, which `git add -A` stages on the next line.
* `git remote remove origin` deletes the tracking refs, then fails on the empty
  `.git/refs/remotes/origin/` directory and prompts. Answering `n` lets it **finish removing
  the config section**, so `git push` then fails with *'origin' does not appear to be a git
  repository*. Inspecting `git remote -v` while that prompt is open is misleading - it still
  shows the remote. Simply `git remote add origin <url>` again; the section is cleanly absent,
  not half-written. The empty directory remains permanently and affects nothing.

---

## What changes

| | Before | After |
|---|---|---|
| `.git` total | 7.5 GB | ~0.15 GB |
| `.git/lfs` | 4.2 GB, 1,768 objects | **19 MB, 19 objects** |
| Tracked files | 806 in the working tree, incl. 479 `.fastq.gz` | **170 files, 134 MB** (158 / 79 MB with the optional trim) |
| Tracked FASTQ | 479 files, 949 MB | **none** - present on disk, ignored |
| GitHub LFS quota | ~1.5 GB of sequence data | ~19 MB of `.rds` / `.RData` |
| Commit history | 3+ commits | one initial commit |

Verified by simulating the fresh repository with `git check-ignore --no-index` against the
current `.gitignore`.

---

## Step 1 - Pause OneDrive

Tray icon -> *Pause syncing* -> 2 hours. OneDrive holds locks on `.git` objects and will
otherwise fight every step below. Do not skip this.

## Step 2 - Set the old history aside (do not delete it yet)

```bat
cd /d "C:\Users\uax75\OneDrive\Documents\R\AAC-BIOcean5D-Coral_Microbiome"
attrib -h .git
move .git _git-old
dir /a
git rev-parse --show-toplevel
```

`move` should report `1 dir(s) moved`, and `git rev-parse` must print
`fatal: not a git repository`. That last line, not the directory listing, is the
authoritative confirmation.

**The `attrib -h` is mandatory in cmd.** Git for Windows marks `.git` as Hidden, and cmd's
`move` refuses hidden files - it fails with `The system cannot find the file specified.`
even though the directory is there (verified by experiment: hidden -> that exact error;
attribute cleared -> `1 dir(s) moved`). For the same reason plain `dir` never lists `.git`,
even though it lists `.gitignore` and `.gitattributes`; use `dir /a`. PowerShell's
`Rename-Item` and bash's `mv` handle hidden directories natively and need no `attrib` step.

Renaming rather than deleting keeps everything up to step 6 reversible: `move _git-old .git`
puts you back exactly where you started. The working tree is untouched either way.

## Step 3 - Initialise a fresh repository, with the git dir OUTSIDE OneDrive

```bat
git lfs version
mkdir C:\git
git init --separate-git-dir="C:\git\AAC-BIOcean5D-Coral_Microbiome.git" -b main
git lfs install --local
```

`--separate-git-dir` keeps the repository where it is - nested with the other R projects and
still backed up by OneDrive - while moving git's object store onto plain local disk. `.git`
becomes a ~150-byte text file containing `gitdir: C:/git/...`; `objects/`, `refs/`, `index`,
`logs/` and `lfs/` all live at the external path. OneDrive then has a tiny text file to sync
instead of several GB of loose objects.

This matters because every obstacle in this repository has come from OneDrive managing
`.git`: loose objects unreadable as cloud placeholders, the `ReadOnly` attribute blocking
`move` with *Access is denied*, reparse tags, and sync locks. OneDrive itself cannot exclude
a folder while keeping it on disk - "Choose folders" selective sync *removes* the local copy,
and "Free up space" only controls hydration - so relocating the git dir is the way to get the
effect of an exclusion.

Two consequences to accept:

* The history is no longer cloud-backed. GitHub is the backup, which is the correct place
  for it.
* The `.git` file holds an **absolute** path. If `C:\git\AAC-BIOcean5D-Coral_Microbiome.git`
  is ever moved or deleted, the working tree stops being a repository until the path in
  `.git` is corrected.

To use an ordinary in-place `.git` instead, substitute `git init -b main`; nothing else in
this runbook changes.

`git lfs version` must report `git-lfs/3.7.0 (GitHub; windows amd64)`. If it reports
something else you are in a shell with a different git, and the LFS checks in step 5 will
not reflect reality.

`.gitattributes` and `.gitignore` are working-tree files, so the new repository picks them
up immediately - the ignore rules are in force for the very first commit.

## Step 4 - Verify BEFORE committing

This is the step that matters.

```bat
git status --porcelain -uall | find /c /v ""
git status --porcelain -uall | findstr _git-old
git status --porcelain -uall | findstr CleanBactSeqs
git status --porcelain -uall | findstr DADAFilteredSeqs
git status --porcelain -uall | findstr fastq.gz
```

Expected: the first prints **158**; the other four print nothing. (It printed 162 on
26 October; the working tree has since gained the section 3.2 window guard and the
re-trimmed environmental inputs. Re-simulated 27 October.)

> **`-uall` is essential.** Without it `git status --porcelain` collapses each untracked
> directory to a single line - it reports 13 lines, not 162 - so the count looks wrong and,
> worse, hides what is inside those directories. `-uall` lists files individually.

> **`_git-old/` must be ignored.** It is a complete 7.5 GB repository sitting inside the
> working tree, so `git add -A` will try to commit all of it: on the first attempt this
> produced a commit of 3,458 files / 7,629 MB, of which 3,296 files were `_git-old`
> (including the 2,872 MB packfile). The rule is now at the top of `.gitignore`. The
> alternative is to move `_git-old` outside the repository instead of ignoring it - but note
> that moving it out of the OneDrive tree forces OneDrive to hydrate any tiered objects
> first, so ignoring it in place and deleting it at step 10 is usually quicker.

Then confirm the large folders are actively ignored - each should echo the matching rule:

```bat
git check-ignore -v "Data/00-AAC-Sequences/CleanBactSeqs/1-TGC01F00POC_S78_L001.fwd_host_removed_R1_001.fastq.gz"
git check-ignore -v "Data/01-AAC-AssignSeqVariants/SilvaDB/silva_nr99_v138.1_train_set.fa.gz"
```

If any FASTQ appears, stop and fix `.gitignore` before continuing.

### The ~55 MB trim (already applied)

Three of the largest remaining files are regenerable intermediates. These lines are already
in `.gitignore` as of 27 October; it takes the first commit from
170 files / 134 MB to 158 files / 79 MB, and LFS from 57 MB to 19 MB:

```
Data/01-AAC-AssignSeqVariants/DADA2-Results/01-AAC-AssignSeqVariants_dadaF.ASV.rds
Data/03-AAC-PhylogeneticTree/FastTree.log
maaslin2_output_directory/
```

To put any of them back in history, delete its line from `.gitignore`.

## Recovery - if you already committed `_git-old` by mistake

Harmless as long as nothing has been pushed. Discard the repository and start again; the
working tree is not touched.

**Rename it rather than deleting it.** The `Everyone:(CI)(DENY)(DC)` ACE on `.git` denies
*delete-child*, which is the right `rmdir /s` needs for entries inside the directory. In
practice `rmdir /s /q .git` clears the deep objects, then stalls near the top with
`The directory is not empty.` followed by `Access is denied.` on the last few entries.
Renaming needs only `DELETE` on the object itself, which you do have:

```bat
cd /d "C:\Users\uax75\OneDrive\Documents\R\AAC-BIOcean5D-Coral_Microbiome"
attrib -h -r .git
move .git _git-dead
```

`_git-dead/` is already in `.gitignore`, so it will not pollute the next commit, and it can
be removed later from a shell with the rights to do so (MobaXterm `rm -rf _git-dead` is the
simplest). Then return to step 3.

Re-running `git add -A` without clearing `.git` would leave several GB of unreferenced
objects behind, so starting from a clean git dir is the right move either way.

If you would rather delete outright, do it from your own elevated shell - not from a tool
session, whose access grants add their own ACEs under `.git`:

```bat
takeown /f .git /r /d y
icacls .git /remove:d Everyone /t
rmdir /s /q .git
```

## Rollback - return to the original repository

At any point before step 6 (deleting the GitHub repository) the whole exercise is
reversible, because the original repository is sitting in `_git-old` untouched.

```bat
cd /d "C:\Users\uax75\OneDrive\Documents\R\AAC-BIOcean5D-Coral_Microbiome"
rem if a new .git exists, move it aside first - it cannot be deleted in place
attrib -h -r .git
move .git _git-dead
rem restore the original
move _git-old .git
git log --oneline -n 3
```

`git log` should show `3372f48`, `96cc524`, `87076e6`. The working tree is never modified by
any of this, so nothing committed or uncommitted is at risk.

`_git-dead` is a husk of empty directories and three 0-4 byte files. It is listed in
`.gitignore`, so it is harmless; remove it whenever convenient from a shell with the rights
(`rm -rf _git-dead` in MobaXterm).

### Run the filesystem steps yourself

Do the renames, deletes and `git init` from your own elevated cmd, not from a tool session.
Agent host-access grants place inheritable ACEs under `.git`, and when that directory is
renamed those ACEs travel with the subtree - which produced repeated *Access is denied*
failures on operations that succeeded perfectly well from an ordinary shell.

---

## Alternative route: orphan commit (no filesystem mutation)

`move .git _git-old` failed with *Access is denied* on both 26 and 27 October, from an
elevated cmd. Diagnosis on 27 October, read-only:

* A freshly created plain directory in the repo root renames fine, so the inherited
  `Everyone:(I)(CI)(DENY)(DC)` on the root is not blocking renames generally.
* `.git` is the only folder carrying that DENY delete-child ACE **explicitly** rather than
  by inheritance, and its OneDrive reparse tag is `0x9000a01a` against `0x9000e01a` on
  ordinary folders - a different cloud-hydration state.
* DENY delete-child explains why `rmdir /s` stalls inside `.git`, but it does not forbid
  renaming `.git` itself. The remaining candidates are the placeholder state and an open
  handle somewhere under `.git`; host processes cannot be enumerated from a tool session, so
  the handle cannot be excluded.

Renaming `.git` was only ever a means to discard history. An orphan commit achieves the same
end without renaming, deleting or moving anything on disk, which sidesteps the blocker
entirely. Verified preconditions on 27 October: refs are `refs/heads/main` plus the two
`refs/remotes/origin/*` tracking refs, **no tags**, and the index holds 618 files of which
479 are `.fastq.gz`.

```bat
cd /d "C:\Users\uax75\OneDrive\Documents\R\AAC-BIOcean5D-Coral_Microbiome"
git checkout --orphan rebuild
git rm -r --cached .
git add -A
git status --porcelain -uall | find /c /v ""
```

`--orphan` starts a branch with no parent but leaves the index and working tree alone, so
`git rm -r --cached .` is required to clear the 618 stale index entries before `git add -A`
re-adds under the current `.gitignore`. Nothing leaves the working tree. The count must be
**158**.

```bat
git commit -m "Rebuild repository without tracked sequence data; add ITS2 (06) and RNAseq (07) analyses and the step-01 amplicon-window guard."
git ls-files | find /c /v ""
git ls-files | findstr fastq.gz
git lfs ls-files | find /c /v ""
```

Expected **158**, nothing, **19**. Then retire the old branch and take its name:

```bat
git branch -D main
git branch -m main
```

Delete and recreate the GitHub repository (runbook steps 6 and 7), then drop the stale
remote so push negotiation cannot be misled by tracking refs that point at commits the
recreated remote has never seen:

```bat
git remote remove origin
git remote add origin https://github.com/EJArmstrong/AAC-BIOcean5D-Coral_Microbiome.git
git push -u origin main
```

### What this leaves behind

GitHub holds one commit and ~19 MB of LFS - the actual objective. The local `.git` stays at
7.5 GB, because the old objects are merely unreferenced, not gone. Reclaiming that space is
cosmetic and can be attempted whenever:

```bat
git reflog expire --expire=now --all
git gc --prune=now
git lfs prune
```

These may fail part-way on the same DENY delete-child ACE. That is harmless - it costs disk,
not correctness. To clear the obstruction first, from your own elevated cmd:

```bat
icacls .git /remove:d Everyone /t
attrib +p .git /s /d
```

`/remove:d Everyone` drops the DENY ACE; `attrib +p` pins the subtree so OneDrive hydrates it
rather than leaving placeholders. Those two lines are also worth one attempt *before* falling
back to this route - if they succeed, `move .git _git-old` may simply work and the main
runbook proceeds unchanged. Close RStudio and any Explorer window on the repository first,
and quit OneDrive rather than merely pausing it.

## Step 5 - Commit

```bat
git add -A
git commit -m "Rebuild repository without tracked sequence data; add ITS2 (06) and RNAseq (07) analyses."
git lfs ls-files | find /c /v ""
git lfs ls-files | findstr fastq.gz
```

Expected: **19** LFS files, all `.rds` / `.RData`, and nothing from the second command.

(Not 20. `.gitattributes` routes `*RData`, a glob that matches files *ending* in `RData`;
`Data/06-AAC-ITS2Symbionts/RData.xlsx` only contains the string and ends in `.xlsx`, so
`git check-attr filter` reports `unspecified` and it is committed as an ordinary blob.
Verified against `git lfs ls-files` on 27 October.)

## Step 6 - Delete the GitHub repository

In the browser: `github.com/EJArmstrong/AAC-BIOcean5D-Coral_Microbiome` -> **Settings** ->
scroll to **Danger Zone** -> *Delete this repository* -> type the name to confirm.

This is the step that actually reclaims the LFS storage, and it is irreversible. It destroys
issues, pull requests, stars, forks and the old commit history. You have confirmed none of
those matter here.

## Step 7 - Recreate it, empty

Still in the browser: **New repository**, same name `AAC-BIOcean5D-Coral_Microbiome`, same
owner, **no** README, **no** `.gitignore`, **no** licence. It must be genuinely empty or the
push will be rejected.

## Step 8 - Push

```bat
git remote add origin https://github.com/EJArmstrong/AAC-BIOcean5D-Coral_Microbiome.git
git push -u origin main
```

## Step 9 - Verify the storage was actually reclaimed

* The repository page should show one commit and no `Data/00-AAC-Sequences/CleanBactSeqs`.
* **github.com -> Settings -> Billing and plans -> Storage -> Git LFS** should fall to
  roughly 19 MB. Check this explicitly; do not assume.

## Step 10 - Clean up locally and resume

Only after the push is verified on GitHub:

```bat
rmdir /s /q _git-old
git status --porcelain
dir /b /s "Data\00-AAC-Sequences\CleanBactSeqs" | find /c /v ""
```

Expected: `git status` prints nothing, and the count is **322** - the sequence data
untouched on disk and now ignored. Then resume OneDrive syncing.

---

## The three numbers to watch

| Checkpoint | Command | Expected |
|---|---|---|
| Before the first commit | `git status --porcelain -uall \| find /c /v ""` | **158** |
| After the commit | `git lfs ls-files \| find /c /v ""` | **19** |
| At the end | `dir /b /s "Data\00-AAC-Sequences\CleanBactSeqs" \| find /c /v ""` | **322** |

## Shell equivalents

Every `git ...` line is identical in any shell; only file operations and the counting
one-liners differ.

| Purpose | Windows cmd | PowerShell | bash / MobaXterm |
|---|---|---|---|
| Enter the repo | `cd /d "<path>"` | `cd '<path>'` | `cd /drives/c/Users/uax75/OneDrive/Documents/R/AAC-BIOcean5D-Coral_Microbiome` |
| Rename `.git` | `attrib -h .git` then `move .git _git-old` | `Rename-Item .git _git-old` | `mv .git _git-old` |
| Undo that rename | `move _git-old .git` | `Rename-Item _git-old .git` | `mv _git-old .git` |
| Show hidden entries | `dir /a` | `Get-ChildItem -Force` | `ls -a` |
| Count lines of output | `... \| find /c /v ""` | `... \| Measure-Object \| %{$_.Count}` | `... \| wc -l` |
| Search output | `... \| findstr PATTERN` | `... \| Select-String 'PATTERN'` | `... \| grep PATTERN` |
| Delete `_git-old` | `rmdir /s /q _git-old` | `Remove-Item -Recurse -Force _git-old` | `rm -rf _git-old` |
| Count files in a folder | `dir /b /s "<dir>" \| find /c /v ""` | `(Get-ChildItem -Recurse -File '<dir>').Count` | `find '<dir>' -type f \| wc -l` |

## Keeping it clean

* `.gitignore` blocks `*.fastq.gz`, `*.fq.gz`, BAM/SAM/CRAM, `.fna`, GFF/GTF, salmon
  indices, `SilvaDB/` and `DADAFilteredSeqs/`. Re-adding any of them needs a deliberate
  `git add -f`.
* `.gitattributes` still routes `*.gz`, `*.rds` and `*RData` to LFS. Keep it: at ~57 MB it
  is doing exactly what LFS is for, and it is the second guard if a `.gz` ever slips past
  `.gitignore`.
* The 41 GB of RNAseq reads and the 347 Mb PMEA_v1 reference stay outside the repository
  permanently - see `Code/07-AAC-RNAseq_quantification.sh`.

## Note

The old commit SHAs (`e30ee61`, `2f8ed51`, `87076e6`, `7b519a2`, `6005d44`) cease to exist
after step 6. Nothing in `Manuscript/` depends on them, but any external note citing a
commit will dangle.
