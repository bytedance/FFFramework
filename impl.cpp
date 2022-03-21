#include <windows.h>
#include <stdint.h>
#define FFFRAMEWORK_EXPORTS
#include "FF_AVC.h"
#include "FF_HEVC.h"
const UINT8* FindRange(const UINT8* s, const UINT8* end)
{
	const UINT8* out = FindCode(s, end);
	if (s < out && out < end && !out[-1])
		out--;
	return out;
}

const INT32 FindCodeSize(const UINT8* s, UINT32 size)
{
	if (size < 4 || s[0] != 0 || s[1] != 0)
		return 0;
	if (s[2] == 1)
		return 3;
	if ((s[2] == 0 && s[3] == 1))
		return 4;
	return 0;
}

void hvcc_init(HEVCDecoderConfigurationRecord* hvcc)
{
	memset(hvcc, 0, sizeof(HEVCDecoderConfigurationRecord));
	hvcc->configurationVersion = 1;
	hvcc->lengthSizeMinusOne = 3;
	hvcc->general_profile_compatibility_flags = 0xffffffff;
	hvcc->general_constraint_indicator_flags = 0xffffffffffff;
	hvcc->min_spatial_segmentation_idc = MAX_SPATIAL_SEGMENTATION + 1;
}

void hvcc_close(HEVCDecoderConfigurationRecord* hvcc)
{
	uint8_t i;
	for (i = 0; i < hvcc->numOfArrays; i++)
	{
		hvcc->array[i].numNalus = 0;
		av_freep(&hvcc->array[i].nalUnit);
		av_freep(&hvcc->array[i].nalUnitLength);
	}
	hvcc->numOfArrays = 0;
	av_freep(&hvcc->array);
}

int get_hvcc(struct HEVCDecoderConfigurationRecord* hvcc, uint8_t* data, size_t size)
{
	const uint8_t* nal_start, * nal_end;
	const uint8_t* end = data + size;
	int type;
	int ret = 0;
	nal_start = FindRange(data, end);
	while (true)
	{
		while (nal_start < end && !*(nal_start++));
		if (nal_start == end)
			break;
		GetBitContext gbc;
		nal_end = FindRange(nal_start, end);
		int size = nal_end - nal_start;
		ret = init_get_bits8(&gbc, nal_start, size);
		if (ret < 0)
			return ret;
		type = (nal_start[0] >> 1) & 0x3F;
		switch (type)
		{
		case MS_HEVC_NAL_VPS:
		case MS_HEVC_NAL_SPS:
		case MS_HEVC_NAL_PPS:
		case MS_HEVC_NAL_SEI_PREFIX:
		case MS_HEVC_NAL_SEI_SUFFIX:
			ret = hvcc_add_nal_unit((uint8_t*)nal_start, size, 1, hvcc);
			if (ret < 0)
				return ret;
			break;
		default:
			ret = AVERROR_INVALIDDATA;
			return ret;
		}
		nal_start = nal_end;
	}
	return ret;
}
