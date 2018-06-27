#include <memory>
#include <string>

namespace spVCF {

struct transcode_stats {
    uint64_t N = 0;                // samples in the project VCF
    uint64_t lines = 0;            // VCF lines (excluding header)
    uint64_t sparse_cells = 0;     // total 'cells' in the sparse representation
    uint64_t sparse90_lines = 0;   // lines encoded with <=10% the dense number of cells
    uint64_t sparse99_lines = 0;   // " <=1% "
    // dense_cells = lines*N

    uint64_t squeezed_cells = 0;   // cells whose QC measures were dropped
};

class Transcoder {
public:
    virtual std::string ProcessLine(const std::string& input_line) = 0;
    virtual transcode_stats Stats() = 0;
};
std::unique_ptr<Transcoder> NewEncoder(bool squeeze);
std::unique_ptr<Transcoder> NewDecoder();

}
