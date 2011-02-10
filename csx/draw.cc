/*
 * draw.cc -- Utility for drawing sparse matrices.
 *
 * Copyright (C) 2009-2011, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2009-2011, Kornilios Kourtis
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */
#include "spm.h"

#include "draw.h"

#include <cairomm/context.h>
#include <cairomm/surface.h>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

namespace bll = boost::lambda;

extern "C"{
	#include <unistd.h>
}

using namespace csx;

void Draw(SPM &spm,
          const char *filename,
		  int row_start, int row_end,
		  const int width, const int height, bool symmetric=false)
{
	//Cairo::RefPtr<Cairo::PdfSurface> surface;
	Cairo::RefPtr<Cairo::ImageSurface> surface;
	Cairo::RefPtr<Cairo::Context> cr;
	SPM::PntIter p, p_start, p_end;
	double max_cols, max_rows, x_pos, y_pos;
	boost::function<double (double v, double max_coord, double max_v)> coord_fn;
	boost::function<double (double a)> x_coord, y_coord;

	//surface = Cairo::PdfSurface::create(filename, width, height);
	surface = Cairo::ImageSurface::create(Cairo::FORMAT_ARGB32, width, height);
	cr = Cairo::Context::create(surface);

	// background
	cr->save();
	cr->set_source_rgb(1, 1, 1);
	cr->paint();
	cr->restore();

	// TODO:
	//  If not HORIZNTAL, transform and transofrm back after end
	if (row_end == 0)
		row_end = spm.GetNrRows();
	max_rows = (double)(row_end - row_start);
	max_cols = (double)spm.GetNrCols();

	// Lambdas for calculationg coordinates
	coord_fn = (((bll::_1 - .5)*(bll::_2))/(bll::_3));
	x_coord = bll::bind(coord_fn, bll::_1, (double)width, max_cols);
	y_coord = bll::bind(coord_fn, bll::_1 - row_start, (double)height, max_rows);

	cr->set_line_cap(Cairo::LINE_CAP_ROUND);
	cr->set_line_width(2.0);
	p_start = spm.PointsBegin(row_start);
	p_end = spm.PointsEnd(row_end);
	for (p = p_start; p != p_end; ++p){
		SpmCooElem elem = *p;
		if (elem.pattern == NULL){
			x_pos = x_coord(elem.x);
			y_pos = y_coord(elem.y);
			//std::cout << elem << "=>" << x_pos << " " << y_pos << "(" << max_rows << ")\n";
			cr->move_to(x_pos, y_pos);
			cr->line_to(x_pos, y_pos);
			if (symmetric){
				x_pos = x_coord(elem.y);
				y_pos = y_coord(elem.x);
				cr->move_to(x_pos, y_pos);
				cr->line_to(x_pos, y_pos);
			}
		} else {
			DeltaRLE::Generator *gen;
			cr->save();
			cr->set_source_rgb(.9, .1, .1);
			gen = elem.pattern->generator(static_cast<CooElem>(elem));
			// FIXME add mapping functions (ontop of {x,y}_coord fneeded
			assert(spm.GetType() == elem.pattern->GetType());
			while ( !gen->IsEmpty() ){
				CooElem e = gen->Next();
				cr->move_to(x_coord(e.x), y_coord(e.y));
				// FIXME: add drawing
			}
			cr->restore();
		}
		cr->stroke();
	}
	//cr->show_page();
	surface->write_to_png(filename);
}

extern int optind;
extern char *optarg;

int main(int argc, char *argv[])
{
	char *out_file = (char *)"out.png";
	bool symmetric = false;
	char c;

	while ((c = getopt(argc, argv, "so:")) != -1){
		switch (c) {
			case 's':
			symmetric = true;
			break;

			case 'o':
			out_file = optarg;
			break;

			default:
			fprintf(stderr, "Error parsing arguments: -%c-\n", c);
			exit(1);

		}
	}

	if (argc - optind < 1){
		std::cout << "Usage: " << argv[0] << " [-s] [-o <out_file>] <mmf_file>\n";
		exit(1);
	}

	SPM *Spm;

	Spm = SPM::LoadMMF(argv[optind]);
	Draw(*Spm, out_file, 0, 0, 800, 800, symmetric);

	return 0;
}
