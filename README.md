# Python_Script_for_Statistical_Resampling_and_Modeling_of_Multivariate_Water_Quality_Data

This special-purpose Python script reads a multivariate water quality data set, generates a larger synthetic data set with matching parameter probability distributions and pairwise correlations, and runs PHREEQC speciation simulations on the synthetic data. The script, as included in this repository, addresses a specific issue at a particular site, but can be generalized/modified to handle similar problems or to set up Monte Carlo geochemical models. The objectives of the script and some example results are shown in more detail in a blog post at [link pending].

The script requires the numpy, pandas, and matplotlib libraries. The following input files are also required:
* components.csv – component/analytes to be listed in the script output; this includes the analyte name and its molecular weight.
* surface_species.txt – a simple list of surface-complexed species, as represented by PHREEQC, to be tracked as “adsorbed”.
* water_quality.csv – observational water quality data set; may include blank entries and comments. “Well”, “Date” and columns must be present, along with corresponding column headings for analytes. The actual example monitoring data used to inform the corresponding post in my blog are omitted here.

Other lists and constraints are included within the script, with comments.

I'd appreciate hearing back from you if you find the code useful. Questions or comments are welcome at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
