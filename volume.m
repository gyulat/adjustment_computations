function vol = volume(x, y, z)
    % Delaunay triangulation in XY
    tri = delaunay(x, y);

    % Vertices of triangles
    X1 = x(tri(:,1)); Y1 = y(tri(:,1)); Z1 = z(tri(:,1));
    X2 = x(tri(:,2)); Y2 = y(tri(:,2)); Z2 = z(tri(:,2));
    X3 = x(tri(:,3)); Y3 = y(tri(:,3)); Z3 = z(tri(:,3));

    % Area of each triangle projected onto XY
    A = 0.5 * abs( X1.*(Y2-Y3) + X2.*(Y3-Y1) + X3.*(Y1-Y2) );

    % Signed volume = area * mean(z)
    vol = sum( A .* (Z1 + Z2 + Z3) / 3 );
end

