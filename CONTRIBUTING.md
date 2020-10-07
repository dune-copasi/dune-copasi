---
id: contribute
title: Contribute
sidebar_label: Contribute
---

When contributing to `DuneCopasi`, first discuss the change you wish to make via the
[issue tracker](https://gitlab.dune-project.org/copasi/dune-copasi/-/issues).
Try to do this before anything else. This will greately improve both your and
our experience working together :-). We follow the
[GitLab Workflow](https://about.gitlab.com/blog/2016/10/25/gitlab-workflow-an-overview/)
in our software development cycle and we expet you to be part of it. In short,
contributing single features can be summarized by following these steps:

1. **ISSUE**: Create an issue in our
   [Gitlab repository](https://gitlab.dune-project.org/copasi/dune-copasi)
   containing a proposal for your modification. We call this a *task*. If your
   proposal has a major change, create a *meta-task*, which should divide the
   proposal in several tasks. Use a suitable issue template for your proposal.
2. In the created issue, discuss with others the whole idea of your proposal.
   This is important because you may get input from someone that know the
   software on how to do it effectively, and poimting out considerations you may
   have missed. It also syncronizes ideas so that changes are more easily
   accepted.
3. **CODE**: Implement what you proposed :-) in a new branch.
4. Create a [Merge Request (MR)](https://docs.gitlab.com/ee/gitlab-basics/add-merge-request.html)
   using suitable template. Mark the MR as
   [*Work In Progres* (WIP)](https://docs.gitlab.com/ee/user/project/merge_requests/work_in_progress_merge_requests.html)
   until is sent to review for approval. MR differ from issues in that MR
   discuss *how* to implement the idea rather than the idea itself. A MR is
   filled of technical details that ensure that the *task* idea is fullfiled.
5. **COMMIT**: In the MR, discuss with others the details of your
   implementation. If discussions show that your proposal has to be modified to
   be accepted, push new commits and update the MR with the new technical
   details.
6. If applicable, create a
   [unit test](https://en.wikipedia.org/wiki/Unit_testing) that proves that your
   implementation has the expected behavior. Ensure that the created test is
   being run by the [Gitlab CI](https://docs.gitlab.com/ee/ci/pipelines.html).
7. **TEST**: Once implementation and test are finished, ensure that your MR is
   passing the different [pipelines](https://gitlab.dune-project.org/copasi/dune-copasi/-/tree/master/.ci).
8. Update the [CHANGELOG](CHANGELOG.md) with changes to the interface. It
   must include a link to the associated MR and a short description that
   summarizes the changes.
9. **REVIEW**: Send the MR to review for approval unmarking the WIP and
   assigning a maintainer to review your modification.

Check out our [code of conduct](CODE_OF_CONDUCT.md), please follow it in our interactions.

