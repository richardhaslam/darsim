function [P,S,z] = Reform_Flash(P,S,z,Grid)

[a,b] = size(P);

if b > 1
    P = reshape(P',Grid.Nx*Grid.Ny,1);
    if (strcmp(S,'Initialize')==1)
    else
        S = reshape(S',Grid.Nx*Grid.Ny,1);
        S(:,2) = 1 - S(:,1);
    end
    if (strcmp(z,'Initialize')==1)
    else
        z = reshape(z',Grid.Nx*Grid.Ny,1);
        z(:,2) = 1 - z(:,1);
    end
else
    P = reshape(P,Grid.Nx,Grid.Ny);
    if (strcmp(S,'Initialize')==1)
    else
        S = reshape(S(:,1),Grid.Ny,Grid.Nx)';
    end
    if (strcmp(z,'Initialize')==1)
    else
        z = reshape(z(:,1),Grid.Ny,Grid.Nx)';
    end
end

end

