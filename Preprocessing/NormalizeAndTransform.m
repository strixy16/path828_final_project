function norm_log_data = NormalizeAndTransform(mat_data)
%NormalizeAndTransform Summary of this function goes here
%   Detailed explanation goes here
    norm_data = SampleNormalizationRF(mat_data);
    norm_log_data = log2(replaceZeros(norm_data, 'lowval'));
end

