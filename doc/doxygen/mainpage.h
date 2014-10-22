/*! \mainpage
 *
 * You will find in this documentation a quick description of mbtrack, a few
 * tutorials, and references.
 *
 * \section intro_sec Introduction
 *
 * \b mbtrack is an electron multi-bunch beam simulation software.
 *
 * It can simulate multiple electrons bunchs circulating in an accelerator's
 * storage ring, and that over many turns.
 *
 * It's purpose is to analyse instabilities growing slowly compared
 * to the revolution period of a bunch in the ring. The step time is defined
 * as this revolution period.
 *
 * \section using Getting started
 *
 * \subsection compile Compile
 *
 * Compile mbtrack using the following commande
 * \code make \endcode
 *
 * \subsection input Create an input
 * Create a input (.conf) file using the <a href="input-example.html">example input file</a>.
 *
 * \subsection run Run
 *
 * Run mbtrack <a href="run-interactively.html">interactively</a> or <a href="run-with-pbs.html">with PBS</a>.
 *
 * \section files Files and directories
 * - \b doc Documentation directory
 *   - \b doc/html/index.html Documentation main page
 * - \b src Source code directory
 *   - \b src/*.c C source files
 *   - \b src/*.h C header files
 * - \b work Work directory, usually used to write simulation input and output files
 * - \b job.mpi PBS job file, used to run mbtrack
 * - \b Makefile Makefile used to build things
 * - \b mbtrack-mpi mbtrack's executable
 */
