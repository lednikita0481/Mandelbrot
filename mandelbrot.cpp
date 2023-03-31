#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include <immintrin.h>

typedef unsigned char BYTE;
const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHTH = 800;
const float X_MAX = 2.0;
const float Y_MAX = 2.0;
const float MAX_DISTANCE = 10.0;
const float DX = 2*X_MAX/WINDOW_WIDTH;
const float DY = 2*Y_MAX/WINDOW_HEIGHTH;

void Draw_Mandelbrot(sf::Texture* Google_Pixel, float x_mov, float y_mov);
void Draw_Mandelbrot_AVX(sf::Texture* Google_Pixel, float x_mov, float y_mov);

int main()
{
    //printf("%d %d %f %f %f %f %f", WINDOW_WIDTH, WINDOW_HEIGHTH, X_MAX, Y_MAX, MAX_DISTANCE, DX, DY);
    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHTH), "Mandelbrot!", sf::Style::Close);
    printf("%d %d %f %f %f %f %f\n", WINDOW_WIDTH, WINDOW_HEIGHTH, X_MAX, Y_MAX, MAX_DISTANCE, DX, DY);
    

    sf::Texture texture;
    texture.create(WINDOW_WIDTH, WINDOW_HEIGHTH);
    Draw_Mandelbrot_AVX(&texture, 0, 0);
    printf("Kuku %d\n", sizeof(sf::Uint8));
    sf::RectangleShape mandelbrot(sf::Vector2f(WINDOW_WIDTH, WINDOW_HEIGHTH));
    mandelbrot.setTexture(&texture);

    float x_mov = 0;
    float y_mov = 0;
    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
            
            float x_mov_old = x_mov;
            float y_mov_old = y_mov;
            
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left))
                x_mov -= 0.1;
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right))
                x_mov += 0.1;
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Up))
                y_mov += 0.1;
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Down))
                y_mov -= 0.1;

            if (x_mov_old != x_mov || y_mov_old != y_mov)
            {   
                if (x_mov >= 2*X_MAX) x_mov -= 2*X_MAX;
                if (y_mov >= 2*Y_MAX) y_mov -= 2*Y_MAX;
                if (x_mov <= -2*X_MAX) x_mov += 2*X_MAX;
                if (y_mov <= -2*Y_MAX) y_mov += 2*Y_MAX;
                Draw_Mandelbrot_AVX(&texture, x_mov, y_mov);
            }


            window.draw(mandelbrot);
            window.display();
        }
    }
}

void Draw_Mandelbrot(sf::Texture* Google_Pixel, float x_mov, float y_mov)
{
    static sf::Uint8 pixels[4*WINDOW_WIDTH*WINDOW_HEIGHTH] = {};
    sf::Image image = Google_Pixel->copyToImage();
    printf("%d %d\n", image.getSize().x, image.getSize().y);

    for (int y0_window = 0; y0_window < WINDOW_HEIGHTH; y0_window++)
    {
        float y0 = (((2*Y_MAX)) * ((float)y0_window/(float)WINDOW_HEIGHTH)) - Y_MAX + y_mov;
        if (y0 > Y_MAX) y0 -= 2*Y_MAX;
        if (y0 < -Y_MAX) y0 += 2*Y_MAX;

        for (int x0_window = 0; x0_window < WINDOW_WIDTH; x0_window++)
        {
            float x0 = (((2*X_MAX))* ((float)x0_window / (float)WINDOW_WIDTH)) - X_MAX - x_mov;

            if (x0 > X_MAX) x0 -= 2*X_MAX;
            if (x0 < -X_MAX) x0 += 2*X_MAX;

            float x = x0;
            float y = y0;
            
            int cur_iter = 1;
            for (; cur_iter < 255; cur_iter++)
            {
                float x_next = (x*x) - (y*y) + x0;
                float y_next = (2*x*y) + y0;

                if ((x_next*x_next) + (y_next*y_next) > (MAX_DISTANCE*MAX_DISTANCE)) break;

                x = x_next;
                y = y_next;
            }

            int i = 4*((WINDOW_WIDTH * y0_window) + x0_window);
            pixels[i] = (sf::Uint8) 2*cur_iter;
            pixels[i+1] = (sf::Uint8) 2*cur_iter;
            pixels[i+2] = (sf::Uint8) 2*cur_iter;
            pixels[i+3] = (sf::Uint8) 255;
        }
    }

    printf("Exit\n");

    Google_Pixel->update(pixels);
    //Google_Pixel->loadFromImage(image);
}

void Draw_Mandelbrot_AVX(sf::Texture* Google_Pixel, float x_mov, float y_mov)
{
    static sf::Uint8 pixels[4*WINDOW_WIDTH*WINDOW_HEIGHTH] = {};


    for (int y0_window = 0; y0_window < WINDOW_HEIGHTH; y0_window++)
    {
        __m256 y0_avx = _mm256_set1_ps((((2*Y_MAX)) * ((float)y0_window/(float)WINDOW_HEIGHTH)) - Y_MAX + y_mov);
        __m256 y_max_avx = _mm256_set1_ps(Y_MAX);
        __m256 min_y_max_avx = _mm256_set1_ps(-Y_MAX);

        __m256 cmp_y0_greater_avx = _mm256_cmp_ps(y0_avx, y_max_avx, _CMP_GT_OQ);
        int cmp_mask_y0_greater = _mm256_movemask_ps(cmp_y0_greater_avx);
        if (cmp_mask_y0_greater != 0)
        {
            y0_avx = _mm256_sub_ps(y0_avx, _mm256_set1_ps(2*Y_MAX));
        }

        __m256 cmp_y0_lower_avx = _mm256_cmp_ps(y0_avx, min_y_max_avx, _CMP_LT_OQ);
        int cmp_mask_y0_lower = _mm256_movemask_ps(cmp_y0_lower_avx);
        if (cmp_mask_y0_lower != 0)
        {
            y0_avx = _mm256_add_ps(y0_avx, _mm256_set1_ps(2*Y_MAX));
        }


        for (int x0_window = 0; x0_window < WINDOW_WIDTH; x0_window+=8)
        {
            __m256 x0_avx = _mm256_set1_ps((((2*X_MAX))* ((float)x0_window / (float)WINDOW_WIDTH)) - X_MAX - x_mov);
            //__m256 x0_avx_offset = _mm256_setr_ps(7*DX, 6*DX, 5*DX, 4*DX, 3*DX, 2*DX, DX, 0);
            __m256 x0_avx_offset = _mm256_setr_ps(0, DX, 2*DX, 3*DX, 4*DX, 5*DX, 6*DX, 7*DX);

            x0_avx = _mm256_add_ps(x0_avx, x0_avx_offset);

            __m256 x_max_avx = _mm256_set1_ps(X_MAX);
            __m256 min_x_max_avx = _mm256_set1_ps(-X_MAX);

            __m256 x_max_avx_2 = _mm256_set1_ps(2*X_MAX);

            __m256 cmp_x0_greater_avx = _mm256_cmp_ps(x0_avx, x_max_avx, _CMP_GT_OQ);
            int cmp_mask_x0_greater = _mm256_movemask_ps(cmp_x0_greater_avx);
            if (cmp_mask_x0_greater != 0)
            {
                x0_avx = _mm256_sub_ps(x0_avx, x_max_avx_2);
            }

            __m256 cmp_x0_lower_avx = _mm256_cmp_ps(x0_avx, min_x_max_avx, _CMP_LT_OQ);
            int cmp_mask_x0_lower = _mm256_movemask_ps(cmp_x0_lower_avx);
            if (cmp_mask_x0_lower != 0)
            {
                x0_avx = _mm256_add_ps(x0_avx, x_max_avx_2);
            }

            __m256 x_avx = x0_avx;
            __m256 y_avx = y0_avx;

            int cur_iter = 0;
            __m256i iterations_avx = _mm256_set1_epi32(0);

            for (; cur_iter < 255; cur_iter++)
            {
                __m256 x_next_avx = _mm256_sub_ps(_mm256_mul_ps(x_avx, x_avx), _mm256_mul_ps(y_avx, y_avx));
                x_next_avx = _mm256_add_ps(x_next_avx, x0_avx);

                __m256 y_next_avx = _mm256_mul_ps(x_avx, y_avx);
                y_next_avx = _mm256_mul_ps(y_next_avx, _mm256_set1_ps(2));
                y_next_avx = _mm256_add_ps(y_next_avx, y0_avx);

                __m256 cmp_iter = _mm256_cmp_ps(_mm256_add_ps(_mm256_mul_ps(x_next_avx, x_next_avx), _mm256_mul_ps(y_next_avx, y_next_avx)), 
                                                _mm256_set1_ps(MAX_DISTANCE*MAX_DISTANCE), _CMP_LT_OS);
                
                int cmp_mask_iter = _mm256_movemask_ps(cmp_iter);

                if (cmp_mask_iter == 0) break;

                __m256i iter_sub = _mm256_castps_si256(cmp_iter);  // 11111111 -> -1 from mask to int (not unsigned, so we sub)

                iterations_avx = _mm256_sub_epi32 (iterations_avx, iter_sub);

                x_avx = x_next_avx;
                y_avx = y_next_avx;
            }

            int* iter_array = (int*) &iterations_avx;
            for (int j = 0; j < 8; j++)
            {
                int i = 4*((WINDOW_WIDTH * y0_window) + x0_window) + 4*j;
                pixels[i] = (sf::Uint8) 2*iter_array[j];
                pixels[i+1] = (sf::Uint8) 2*iter_array[j];
                pixels[i+2] = (sf::Uint8) 2*iter_array[j];
                pixels[i+3] = (sf::Uint8) 255;
            }

        }

    }

    Google_Pixel->update(pixels);
}