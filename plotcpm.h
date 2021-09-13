#ifndef PLOTCPM_H
#define PLOTCPM_H
#include <QString>
#include "qtgraph.h"
#include "structures.h"
#include "functions.h"
#include "myparameters.h"

class PlotCPM : public QtGraphics
{
    Q_OBJECT
public:
    explicit PlotCPM(int xfield, int yfield, char *moviefile=0);
    void Plot(VOX* pv);
    void PlotHueStrainField(bool STRAINFIELD);
    void PlotPrincipleStrainField(VOX *pv);
    double MaxStrainMagnitude(void) const;
    double MaxForceMagnitude(void) const;
    void CalculateStrainField(NOD *pn) const;
    void CalculateForceField(NOD *pn) const;
    void PlotNodalForces(NOD *pn);
    void PlotNodeConnection(NOD *pn, int *pathx, int *pathy);
    void StrainColorBar(void);
    void DrawColorBarLabel(int yval, double val);
signals:
    
public slots:
    
private:
    double ***strain;
    double ***force;
};

#endif // PLOTCPM_H
