<p align="center">
<img src="../media/outline.png" width="400">
</p>


# Introduction

As [explained by George Whitesides](../media/whiteside.pdf), an _outline_ is a system for planning your in-progress research. It consists of carefully orga­nized and presented data, with attendant objectives, hypotheses, and conclusions.

An outline itself contains little text as it will be rewritten many times throughout the course of your research as part of our continuous effort to un­derstand, analyze, summarize, and reformulate hypotheses on paper. In order to manage our research projects, we augment each outline with a ZenHub project board that gives an overview of current and planned issues in a [Kanban-like framework](https://www.atlassian.com/agile/kanban).


# The outline repository template

* [`code`](code/README.md) - All code that is required to reproduce the results of this project.
* [`data`](data/README.md) - All input/output data that is relevant to this project.
* [`outline`](outline/README.md) - The actual <G|QC|J> outline.
* [`notes`](notes/README.md) - A set of notes that explain context and background of this research notes.


# DevOps

## Devcontainers / Docker

The `outline` folder contains its own `.devcontainer`, such that you should be able to compile the outline on your infrastructure by opening **only** the `outline` folder in VS Code. (When opening the `outline` folder in VS Code, it will prompt you to re-open the folder in that container.) Similarly, the `notes` folder also contains its own `.devcontainer`.

We urge you to design a separate `.devcontainer` for the source code you've written, such that all members involved in your research can reproduce your results on their infrastructure.


## GitHub Actions

This repo is under [CI/CD with Github Actions](https://github.com/features/actions). As such, each contributing member should be able to download the latest version of your outline and notes by clicking on the green :heavy_check_mark: next to that commit > Details > Artifacts (top right) > Download artifacts.

# Publishing: when outlines become research papers

When your outline has matured through several iterations, you can start communicating your research to the general public outside <G|QC|G> by trying to [publish your research](https://publish.acs.org/publish/author_university), which is an [art](https://pubs.acs.org/page/vi/art_of_scientific_publication.html) in itself.

## Journals of interest

### Main journals

- [PCCP](https://www.rsc.org/journals-books-databases/about-journals/pccp/) - focused on important results in the area of physical chemistry. When to choose this journal?: key chemical insights based on analyses from quantum chemical data.
- [JCTC](https://pubs.acs.org/journal/jctcce) - focused new theories, methodology, and/or important applications in quantum electronic structure, molecular dynamics, and statistical mechanics. When to choose this journal?: novel quantum chemical theories / algorithms aimed at the developer community.
- [JCC](https://onlinelibrary.wiley.com/journal/1096987x) - concerned with all aspects of computational chemistry. When to choose this journal?: computational studies that apply existing techniques in new context.
- [JPC-A](https://pubs.acs.org/journal/jpcafh) - focuses on novel tools and methods in experiment and theory. When to choose this journal?: new tools or methods that are of broad interest to the physical chemistry community.
- [JCP](https://aip.scitation.org/journal/jcp) - focuses on advanced theoretical and computational methods of relevance to essentially all branches of chemical physics. Papers on software packages that implement chemical physics methods are also welcome. When to choose this journal?: new tools or methods that are of broad interest to the physical chemistry community.

### Top-tier journals

- [WIREs Computational Molecular Science](https://onlinelibrary.wiley.com/journal/17590884) -  promotes cross-disciplinary research on computational chemistry, biochemistry and materials science and provides authoritative, encyclopedic resources addressing key topics from diverse research perspectives. When to choose this journal?: new insights obtained by reframing and elucidating the links between a large number of works, some of which originated from <G|QC|G>.
- [Nature Chemistry](https://www.nature.com/nchem/) - focuses on publishing high-quality papers that describe the most significant and cutting-edge research in all areas of chemistry. When to choose this journal?: ground breaking work that is expected to have a large impact in chemistry.
- [PRL](https://journals.aps.org/prl/about) - focuses on short, high-quality reports of the most influential developments and transformative ideas in the full arc of fundamental and interdisciplinary physics research. When to choose this journal?: ground breaking work that is expected to have a large impact in chemistry and physics.

## Cover letters

- [Writing cover letters for scientific manuscripts](https://biosciencewriters.com/Writing-Cover-Letters-for-Scientific-Manuscripts.aspx)
- [Nature: how to write a cover letter](http://blogs.nature.com/methagora/2013/09/how-to-write-a-cover-letter.html)
- [Nature immunology: prelude to a good story](https://www.nature.com/articles/ni0208-107)
- [Springer: cover letters](https://www.springer.com/gp/authors-editors/authorandreviewertutorials/submitting-to-a-journal-and-peer-review/cover-letters/10285574)
- [How to write the best journal submission cover letter](https://wordvice.com/journal-submission-cover-letter/)

## Rebuttal letters

- [How to write a rebuttal letter](http://blogs.nature.com/methagora/2013/09/how-to-write-a-rebuttal-letter.html)
