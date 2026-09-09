# Cross-platform installation evaluation (2026-09-08)

## Scope

This evaluation retested the supported Windows portable-Cygwin and Ubuntu
installation paths, the installer smoke test, focused regression tests, and a
real PGC schizophrenia sex-stratified local-locus workflow. The real-data test
used `rs185665940,rs10166057,rs4852780`, with `rs10166057` as the PLINK2 LD
reference and an R2 threshold of zero.

## Environments

- Windows 10 portable Cygwin, Perl architecture
  `x86_64-cygwin-threads-multi`, gnuplot 6.0 patchlevel 5.
- Ubuntu 24.04.4 LTS clean checkout in an isolated WSL1 distribution, Perl
  architecture `x86_64-linux-gnu-thread-multi`, gnuplot 6.0 patchlevel 0.
  The Canonical minimal rootfs SHA256 was verified before import.
- Vagrant 2.4.9: `vagrant validate` passed for
  `install/vagrant/ubuntu/Vagrantfile`.

The workstation could not boot the Hyper-V Vagrant VM from the available
non-administrative shell, and firmware virtualization was unavailable to WSL2
and Docker Desktop. Consequently, Ubuntu execution was validated in a clean
Ubuntu 24.04 userland instead of claiming that a Vagrant VM boot succeeded.

## Defects found and corrected

1. Passing 65,405 SNP and R2 values as command-line arguments exceeded Linux
   `ARG_MAX`. The wrapper now writes a per-locus `*.ld_r2.tsv` sidecar and the
   renderer accepts `--ld-r2-file`.
2. `DBI`, `DBD::SQLite`, and `Text::CSV` were used by the SQLite HaploReg cache
   tools but absent from `cpanfile`. They are now installed and checked.
3. The two SQLite-backed helper scripts did not discover repository-local Perl
   modules in a fresh shell. They now bootstrap the platform-specific local
   Perl tree themselves.
4. An in-shell Cygwin package refresh requested a global runtime upgrade,
   causing locked-DLL failures. The in-shell installer now adds required
   packages without `-g`; the external PowerShell bootstrap remains responsible
   for safe full refreshes after stopping portable-Cygwin processes.
5. The obsolete Cygwin package name `openssl-devel` was replaced with
   `libssl-devel`. Prebuilt `perl-DBI`, `perl-DBD-SQLite`, and `perl-Text-CSV`
   packages are preferred to slow native CPAN builds.
6. Vagrant previously installed Linux native files directly into the Windows
   shared checkout. Provisioning now copies source into
   `/home/vagrant/MultiGWAS-Explorer` and installs on the VM filesystem.
7. Repository attributes now keep shell, Perl, Python, Docker, and Vagrant
   sources LF-terminated while preserving CRLF for Windows launcher files.

## Automated checks

The following passed on both portable Cygwin and Ubuntu:

- `install/check_pipeline_install.sh`
- Perl syntax checks for the main SAS and gnuplot entry points
- pairwise GWAS allele harmonization
- short and long QC report fixtures
- gnuplot directory-layout test
- signed-R2 heatmap rendering and SAS/gnuplot contract tests
- SAS ODA debug-macro and space-exhaustion guards
- SQLite HaploReg cache LD-clumping fixture
- synthetic top-hit MAF safeguards

`vagrant validate` also passed after the VM-local installation change. An
actual SAS ODA login/submission was not used as an Ubuntu installation test;
that operation requires the user's ODA credentials and an available remote SAS
session.

## Real PGC schizophrenia validation

The three-SNP local Manhattan, local GTF, and forest workflow completed on the
Ubuntu installation after the file-based LD handoff. The local-GTF data
contained 37,377 plotted rows; 37,026 had numeric LD R2 values spanning 0 to 1.
The signed color values spanned -1 to +1 (17,726 positive, 18,720 negative,
and 931 zero).

The same file-based LD renderer was also exercised successfully from portable
Cygwin. The generated PGC local-GTF HTML and PNG are under the output directory
configured by `configs/spec_pgc_scz_sex_common_automation.json`.

## Recommended rerun

Windows PowerShell:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass `
  -File .\install\install_windows_portable_cygwin.ps1
```

Ubuntu:

```bash
sudo bash install/install_ubuntu.sh
bash install/check_pipeline_install.sh
```

Vagrant on a host with hardware virtualization and an administrative Hyper-V
session:

```powershell
cd install\vagrant\ubuntu
vagrant up --provider=hyperv
vagrant ssh -c "cd /home/vagrant/MultiGWAS-Explorer && bash install/check_pipeline_install.sh"
```
