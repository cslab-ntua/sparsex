#include "spm.h"

#include "draw.h"

#include <cairomm/context.h>
#include <cairomm/surface.h>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

namespace bll = boost::lambda;

using namespace csx;

void Draw(SPM &spm,
          const char *filename,
		  int row_start, int row_end,
		  const int width, const int height)
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
		row_end = spm.nrows;
	max_rows = (double)(row_end - row_start);
	max_cols = (double)spm.ncols;

	// Lambdas for calculationg coordinates
	coord_fn = (((bll::_1 - .5)*(bll::_2))/(bll::_3));
	x_coord = bll::bind(coord_fn, bll::_1, (double)width, max_cols);
	y_coord = bll::bind(coord_fn, bll::_1 - row_start, (double)height, max_rows);

	p_start = spm.points_begin(row_start);
	p_end = spm.points_end(row_end);
	for (p = p_start; p != p_end; ++p){
		SpmCooElem elem = *p;
		if (elem.pattern == NULL){
			x_pos = x_coord(elem.x);
			y_pos = y_coord(elem.y);
			//std::cout << elem << "=>" << x_pos << " " << y_pos << "(" << max_rows << ")\n";
			cr->move_to(x_pos, y_pos);
			cr->show_text("x");
		} else {
			Pattern::Generator *gen;
			cr->save();
			cr->set_source_rgb(.9, .1, .1);
			gen = elem.pattern->generator(static_cast<CooElem>(elem));
			// FIXME add mapping functions (ontop of {x,y}_coord fneeded
			assert(spm.type == elem.pattern->type);
			while ( !gen->isEmpty() ){
				CooElem e = gen->next();
				cr->move_to(x_coord(e.x), y_coord(e.y));
				cr->show_text("x");
			}
			cr->restore();
		}
		cr->stroke();
	}
	//cr->show_page();
	surface->write_to_png(filename);
}

int main(int argc, const char *argv[])
{
	if (argc < 2){
		std::cout << "Usage: " << argv[0] << " <mmf_file> [<out_file>]\n";
		exit(1);
	}

	const char *out_file = argc > 2 ? argv[2] : "test.png";

	SPM *Spm;

	Spm = SPM::loadMMF(argv[1]);
	Draw(*Spm, out_file, 0, 0, 600, 800);

	return 0;
}
