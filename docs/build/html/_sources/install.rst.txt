Howto Install CASCADe
=====================

### Clone a repository

To start working locally on an existing remote repository,
clone it with the command `git clone <repository path>`.
By cloning a repository, you'll download a copy of its
files into your local computer, preserving the Git
connection with the remote repository.

You can either clone it via HTTPS or [SSH](../ssh/README.md).
If you chose to clone it via HTTPS, you'll have to enter your
credentials every time you pull and push. With SSH, you enter
your credentials once and can pull and push straightaway.

You can find both paths (HTTPS and SSH) by navigating to
your project's landing page and clicking **Clone**. GitLab
will prompt you with both paths, from which you can copy
and paste in your command line.

As an example, consider a repository path:

- HTTPS: `https://gitlab.com/jbouwman/CASCADe`
- SSH: `` git@gitlab.com:jbouwman/CASCADe ``

To get started, open a terminal window in the directory
you wish to clone the repository files into, and run one
of the following commands.

Clone via HTTPS:

```bash
git clone https://gitlab.com/jbouwman/CASCADe
```

Clone via SSH:

```bash
git clone git@gitlab.com:jbouwman/CASCADe 
```

Both commands will download a copy of the files in a
folder named after the project's name.

You can then navigate to the directory and start working
on it locally.

### Go to the master branch to pull the latest changes from there

```bash
git checkout master
```

### Download the latest changes in the project

This is for you to work on an up-to-date copy (it is important to do this every time you start working on a project), while you set up tracking branches. You pull from remote repositories to get all the changes made by users since the last time you cloned or pulled the project. Later, you can push your local commits to the remote repositories.

```bash
git pull REMOTE NAME-OF-BRANCH
```

When you first clone a repository, REMOTE is typically "origin". This is where the repository came from, and it indicates the SSH or HTTPS URL of the repository on the remote server. NAME-OF-BRANCH is usually "master", but it may be any existing branch.

